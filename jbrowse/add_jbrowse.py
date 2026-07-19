#!/usr/bin/env python3
"""
add_jbrowse.py — register a JBrowse2 genome for the trackplot 'JB' button.

STANDALONE config helper. It has NOTHING to do with the bscore pipeline
(pval_dist.py); it only edits JSON so the app's per-window 'JB' button can open a
genome browser at the selected locus. It touches exactly three files:

  jbrowse/config.json      : add an assembly (remote reference streamed via the
                             proxy) + a gene-model track + any local dataset tracks
  jbrowse/proxies.json     : add  /ref/<name>/ -> <proxy-base>.  serve.py merges this
                             at startup and adds CORS + a browser User-Agent, so hosts
                             that block python-urllib or send no CORS still stream.
  bin/data/datasets.json   : set this dataset's "jbrowse" URL (what the JB button opens)

Remote genomes are streamed by HTTP Range — never downloaded. The host only needs to
serve the .fa (+ .fai) or .fa.gz (+ .fai + .gzi) and honour Range requests.

The command is idempotent: re-running for the same --assembly/--name replaces the
prior assembly, its tracks, and its proxy entry.

Examples
--------
Arabidopsis (TAIR JBrowse2, plain .fa + remote .fai + tabixed gene models):
  python3 add_jbrowse.py --name aras --assembly TAIR10 \
      --proxy-base https://jbrowse2.arabidopsis.org/ \
      --fasta TAIR9_chr_all.fa --fai TAIR9_chr_all.fa.fai \
      --gff q4-2024/GL.gff3.gz --gff-tbi q4-2024/GL.gff3.gz.tbi \
      --aliases TAIR9 Athaliana_447_TAIR10

Bgzip-compressed reference (needs --fai and --gzi):
  python3 add_jbrowse.py --name rh --assembly RH_v3 \
      --proxy-base https://spuddb.uga.edu/jb2/ \
      --fasta RH_v3.asm.fa.gz --fai RH_v3.asm.fa.gz.fai --gzi RH_v3.asm.fa.gz.gzi \
      --gff gff/RH_v3.working_models.jb2.sort.gff3.gz \
      --gff-tbi gff/RH_v3.working_models.jb2.sort.gff3.gz.tbi
"""
import argparse, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "config.json")
PROXIES = os.path.join(HERE, "proxies.json")
DATASETS = os.path.normpath(os.path.join(HERE, "..", "bin", "data", "datasets.json"))


def load(path, default):
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    return default


def save(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
        f.write("\n")


def uri(u):
    return {"locationType": "UriLocation", "uri": u}


def main():
    ap = argparse.ArgumentParser(
        description="Register a JBrowse2 genome for the trackplot 'JB' button (edits JSON only).")
    ap.add_argument("--name", required=True,
                    help="dataset name; matches bin/data/<name>/ and the datasets.json entry")
    ap.add_argument("--assembly", required=True, help="JBrowse assembly name, e.g. TAIR10")
    ap.add_argument("--proxy-base", required=True,
                    help="remote host base URL the paths hang off, e.g. https://jbrowse2.arabidopsis.org/")
    ap.add_argument("--fasta", required=True, help="reference path under the proxy base (.fa or .fa.gz)")
    ap.add_argument("--fai", help="index path under the proxy base; omit to use local tracks/<basename>.fai")
    ap.add_argument("--gzi", help="bgzip .gzi path (REQUIRED when --fasta ends in .gz)")
    ap.add_argument("--gff", help="gene-model GFF3 path under the proxy base (optional)")
    ap.add_argument("--gff-tbi", help="tabix index for --gff (implies Gff3TabixAdapter; else whole-file Gff3Adapter)")
    ap.add_argument("--aliases", nargs="*", default=[], help="assembly aliases")
    ap.add_argument("--port", default="9000", help="serve.py port used in the jbrowse URL (default 9000)")
    ap.add_argument("--data-dir", default=os.path.normpath(os.path.join(HERE, "..", "bin", "data")),
                    help="dataset data dir, for auto-discovering local tracks (default ../bin/data)")
    ap.add_argument("--no-local-tracks", action="store_true",
                    help="do not add transcript-region / slice-site / read tracks from bin/data/<name>/")
    args = ap.parse_args()

    name, asm = args.name, args.assembly
    prefix = f"/ref/{name}/"
    base = args.proxy_base if args.proxy_base.endswith("/") else args.proxy_base + "/"

    # ---- 1) proxy prefix -> proxies.json (serve.py merges it) --------------------
    proxies = load(PROXIES, {})
    proxies[prefix] = base
    save(PROXIES, proxies)

    # ---- 2) assembly + tracks -> config.json ------------------------------------
    cfg = load(CONFIG, {"assemblies": [], "tracks": []})
    cfg.setdefault("assemblies", []); cfg.setdefault("tracks", [])

    if args.fasta.endswith(".gz"):
        if not args.fai or not args.gzi:
            sys.exit("error: a .fa.gz reference needs both --fai and --gzi (BgzipFastaAdapter)")
        adapter = {"type": "BgzipFastaAdapter",
                   "fastaLocation": uri(prefix + args.fasta),
                   "faiLocation":   uri(prefix + args.fai),
                   "gziLocation":   uri(prefix + args.gzi)}
    else:
        fai = (prefix + args.fai) if args.fai else f"tracks/{os.path.basename(args.fasta)}.fai"
        adapter = {"type": "IndexedFastaAdapter",
                   "fastaLocation": uri(prefix + args.fasta),
                   "faiLocation":   uri(fai)}

    cfg["assemblies"] = [a for a in cfg["assemblies"] if a.get("name") != asm]
    cfg["assemblies"].append({
        "name": asm, "aliases": args.aliases,
        "sequence": {"type": "ReferenceSequenceTrack",
                     "trackId": f"{asm}-ReferenceSequenceTrack", "adapter": adapter},
    })

    # rebuild every track pointing at this assembly
    cfg["tracks"] = [t for t in cfg["tracks"] if asm not in t.get("assemblyNames", [])]
    added, prep = [], []

    if args.gff:
        if args.gff_tbi:
            gm = {"type": "Gff3TabixAdapter", "gffGzLocation": uri(prefix + args.gff),
                  "index": {"location": uri(prefix + args.gff_tbi)}}
        else:
            gm = {"type": "Gff3Adapter", "gffLocation": uri(prefix + args.gff)}
        cfg["tracks"].append({
            "type": "FeatureTrack", "trackId": f"{name}_gm", "name": f"{name} gene models",
            "assemblyNames": [asm], "category": ["Gene Annotation"], "adapter": gm,
            "displays": [{"type": "LinearBasicDisplay", "displayId": f"{name}_gm-LinearBasicDisplay"},
                         {"type": "LinearArcDisplay", "displayId": f"{name}_gm-LinearArcDisplay"}],
        })
        added.append(f"{name}_gm")

    if not args.no_local_tracks:
        dsdir = os.path.join(args.data_dir, name)

        def find(sub, needles, avoid=()):
            d = os.path.join(dsdir, sub)
            if not os.path.isdir(d):
                return None
            for fn in sorted(os.listdir(d)):
                low = fn.lower()
                if all(n in low for n in needles) and not any(a in low for a in avoid):
                    return os.path.join(d, fn)
            return None

        tr = find("genome", ["transcript", "region", ".bed"])
        if tr:
            b = os.path.basename(tr)
            cfg["tracks"].append({
                "type": "FeatureTrack", "trackId": f"{name}_transcript_regions",
                "name": f"PPR transcript regions ({name})", "assemblyNames": [asm],
                "category": ["PPR"],
                "adapter": {"type": "BedAdapter", "bedLocation": uri(f"tracks/{b}")}})
            added.append(f"{name}_transcript_regions")
            prep.append(f'cp "{tr}" tracks/{b}')

        sl = find(os.path.join("slice", "slice_bed"), ["ft", ".bed"], avoid=["origin"])
        if sl:
            cfg["tracks"].append({
                "type": "FeatureTrack", "trackId": f"{name}_slice_sites",
                "name": f"Predicted slice sites ({name})", "assemblyNames": [asm],
                "category": ["sRNA", "targets"],
                "adapter": {"type": "BedTabixAdapter",
                            "bedGzLocation": uri(f"tracks/{name}_slice_sites.sorted.bed.gz"),
                            "index": {"indexType": "TBI",
                                      "location": uri(f"tracks/{name}_slice_sites.sorted.bed.gz.tbi")}}})
            added.append(f"{name}_slice_sites")
            prep.append(f'LC_ALL=C sort -k1,1 -k2,2n "{sl}" | bgzip > tracks/{name}_slice_sites.sorted.bed.gz '
                        f'&& tabix -p bed tracks/{name}_slice_sites.sorted.bed.gz')

        bam = find("bam", [".bam"], avoid=[".bai", ".csi"])
        if bam:
            cfg["tracks"].append({
                "type": "FeatureTrack", "trackId": f"{name}_reads",
                "name": f"sRNA alignments ({name})", "assemblyNames": [asm],
                "category": ["sRNA"],
                "adapter": {"type": "BedTabixAdapter",
                            "bedGzLocation": uri(f"tracks/{name}_reads.sorted.bed.gz"),
                            "index": {"indexType": "TBI",
                                      "location": uri(f"tracks/{name}_reads.sorted.bed.gz.tbi")}}})
            added.append(f"{name}_reads")
            prep.append(f'bedtools bamtobed -i "{bam}" | LC_ALL=C sort -k1,1 -k2,2n | bgzip > '
                        f'tracks/{name}_reads.sorted.bed.gz && tabix -p bed tracks/{name}_reads.sorted.bed.gz')

    save(CONFIG, cfg)

    # ---- 3) jbrowse URL -> datasets.json ----------------------------------------
    tracks_param = f"{name}_gm" if args.gff else (added[0] if added else "")
    jburl = f"http://localhost:{args.port}/?assembly={asm}" + (f"&tracks={tracks_param}" if tracks_param else "")
    ds = load(DATASETS, None)
    if ds and "datasets" in ds:
        hit = next((d for d in ds["datasets"] if d.get("name") == name), None)
        if hit is not None:
            hit["jbrowse"] = jburl
            save(DATASETS, ds)
        else:
            print(f"note: no '{name}' entry in datasets.json — set its \"jbrowse\" to:\n  {jburl}")
    else:
        print(f"note: datasets.json not found — set the dataset's \"jbrowse\" to:\n  {jburl}")

    # ---- report -----------------------------------------------------------------
    print(f"registered assembly '{asm}'  (proxy {prefix} -> {base})")
    print(f"  tracks: {', '.join(added) if added else '(none)'}")
    print(f"  JB url: {jburl}")
    if prep:
        print("\nLocal tracks are referenced but not yet built. From jbrowse/, run:")
        print("  mkdir -p tracks")
        for c in prep:
            print("  " + c)
    print("\nRestart serve.py so the new proxy takes effect:  python3 serve.py " + args.port)


if __name__ == "__main__":
    main()
