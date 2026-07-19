#include "ofApp.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>

static const char* kDir =
    "/Users/vef25hok/Downloads/ppr_srna_local/trackplot/testdata/rh/";

static std::vector<std::string> splitTab(const std::string& s) {
    std::vector<std::string> out; std::string cur;
    for (char c : s) { if (c == '\t') { out.push_back(cur); cur.clear(); } else cur.push_back(c); }
    out.push_back(cur); return out;
}

// fast ASCII lowercase (IDs are ASCII) — avoids ofToLower's per-call locale/UTF8 allocation
static std::string lc(const std::string& s) {
    std::string r = s;
    for (char& c : r) if (c >= 'A' && c <= 'Z') c = (char)(c + 32);
    return r;
}

static void drawArrowAt(float ang, float radius, float arrowR) {
    float tipR = radius + arrowR - 16, baseR = radius + arrowR, half = 0.025f;
    ofDrawTriangle(tipR*cosf(ang), tipR*sinf(ang),
                   baseR*cosf(ang-half), baseR*sinf(ang-half),
                   baseR*cosf(ang+half), baseR*sinf(ang+half));
}

void ofApp::setup() {
    ofSetWindowTitle("ppr circos");
    ofSetEscapeQuitsApp(false);
    ofSetCircleResolution(64);
    ofEnableAntiAliasing();

    buildDatasets();

    // GW live phase/DC helper: python with pysam + the gw_phase.py script (dataDir/../../)
    gwPython = ofFile::doesFileExist("/opt/homebrew/bin/python3.14") ? "/opt/homebrew/bin/python3.14" : "python3";
    gwScript = ofFilePath::getAbsolutePath(ofFilePath::join(dataDir, "../../gw_phase.py"));
    if (!ofFile::doesFileExist(gwScript))
        gwScript = "/Users/vef25hok/Downloads/ppr_srna_local/trackplot/gw_phase.py";

    datasetParam.set("dataset", "");
    datasetDropdown = std::make_unique<ofxDropdown>(datasetParam);
    for (auto& d : datasets) datasetDropdown->add(d.name);
    datasetDropdown->disableMultipleSelection();
    gui.setup("dataset");
    gui.setPosition(ofGetWidth() - 220, 12);
    gui.add(datasetDropdown.get());

    datasetParam.addListener(this, &ofApp::onDatasetChanged);
    datasetParam = datasets[0].name;      // fires onDatasetChanged -> switchDataset(0)
}

void ofApp::onDatasetChanged(std::string& name) {
    for (int i = 0; i < (int)datasets.size(); ++i)
        if (datasets[i].name == name) { switchDataset(i); return; }
}

void ofApp::buildDatasets() {
    datasets.clear();
    std::string base = dataDir.empty() ? ofToDataPath("", true) : dataDir;   // default: bin/data
    if (!base.empty() && base.back() != '/') base += '/';
    dataDir = base;
    ofLogNotice("ofApp") << "data dir: " << base;

    // Option B: read datasets from <dataDir>/datasets.json (written by pval_dist.py --install)
    std::string manifest = base + "datasets.json";
    if (ofFile::doesFileExist(manifest)) {
        ofJson j = ofLoadJson(manifest);
        for (auto& d : j["datasets"]) {
            Dataset ds;
            ds.name = d.value("name", std::string());
            if (d.contains("bw")) for (auto& b : d["bw"]) ds.bw.push_back(base + b.get<std::string>());
            ds.tsv         = d.contains("bscore")      ? base + d["bscore"].get<std::string>()      : "";
            ds.origins     = d.contains("origins")     ? base + d["origins"].get<std::string>()     : "";
            ds.fasta       = d.contains("fasta")       ? base + d["fasta"].get<std::string>()       : "";
            ds.transcript  = d.contains("transcript")  ? base + d["transcript"].get<std::string>()  : "";
            ds.fai         = d.contains("fai")         ? base + d["fai"].get<std::string>()         : "";
            ds.clustersOnTx= d.contains("clustersOnTx")? base + d["clustersOnTx"].get<std::string>(): "";
            if (d.contains("bam")) {                                   // absolute path used as-is
                std::string b = d["bam"].get<std::string>();
                ds.bam = (!b.empty() && b[0] == '/') ? b : base + b;
            }
            ds.jbrowse = d.value("jbrowse", std::string());            // URL, not a path -> no base prefix
            datasets.push_back(ds);
        }
        ofLogNotice("ofApp") << "loaded " << datasets.size() << " datasets from " << manifest;
        if (!datasets.empty()) return;
    }

    // fallback (pre-migration): the old hardcoded testdata layout
    ofLogWarning("ofApp") << "no datasets.json — using hardcoded testdata fallback";
    std::string RH  = "/Users/vef25hok/Downloads/ppr_srna_local/trackplot/testdata/rh/";
    std::string DMC = "/Users/vef25hok/Downloads/ppr_srna_local/trackplot/testdata/dm/ctrl_dm_clusflank/";
    std::string DMI = "/Users/vef25hok/Downloads/ppr_srna_local/trackplot/testdata/dm/infec_dm_clusflank/";
    datasets = {
        { "RH baldrich",
          { RH + "rhppr_bw_fwd.bw", RH + "rhppr_bw_rev.bw" },
          RH + "rh_bscore.tsv", RH + "trigger_origins.bed",
          RH + "mir_and_baldrich_uncoll_condensedmin20_21_23.fa", "" },
        { "DM ctrl",
          { DMC + "dm_ctrl_bw/dm_ctrl_bw_fwd.bw", DMC + "dm_ctrl_bw/dm_ctrl_bw_rev.bw" },
          DMC + "dm_ctrl_bscore.tsv", DMC + "trigger_origins.bed",
          DMC + "cleaveland/mir_and_merged_alignmentsmin20_21_23.fa",
          DMC + "dm_pare_pinfes_ctrl_ppregion_clusters_on_tx.bed" },
        { "DM infec",
          { DMI + "dm_infec_bw/dm_infec_bw_fwd.bw", DMI + "dm_infec_bw/dm_infec_bw_rev.bw" },
          DMI + "dm_infec_bscore.tsv", DMI + "trigger_origins.bed",
          DMI + "cleaveland/mir_and_merged_alignmentsmin20_21_23.fa",
          DMI + "dm_pare_pinfes_infec_ppregion_clusters_on_tx.bed" },
    };
}

void ofApp::switchDataset(int i) {
    if (i < 0 || i >= (int)datasets.size()) return;
    curDataset = i;
    const Dataset& d = datasets[i];
    bwPaths = d.bw; tsvPath = d.tsv; originsPath = d.origins; fastaPath = d.fasta;
    jbUrl = d.jbrowse;                                 // JBrowse URL template for the 'JB' button

    // reset view state
    selRow = -1; selOrigins.clear(); tableRows.clear(); activeChrom.clear();
    selectedArc = -1; mode = Mode::Overview; rotation = 0.f; tableScroll = 0;

    openReaders(); loadChroms(); layoutArcs(); loadTrack(300);
    loadTsv(tsvPath); loadOrigins(originsPath); loadFasta(fastaPath);
    loadClusterTx(datasets[i].clustersOnTx);

    // dataset dir = <dataDir>/<name>/ ; auto-discover transcript BED (genome/) and BAM (bam/)
    std::string dsDir = dataDir + d.name + "/";

    // transcript-regions BED in <dataset>/genome/ (datasets.json "transcript" overrides)
    std::string txPath = d.transcript;
    if (txPath.empty()) {                                  // only a transcript/region-named BED (NOT a stray alignments BED)
        ofDirectory gdir(dsDir + "genome"); gdir.allowExt("bed"); gdir.listDir();
        for (std::size_t k = 0; k < gdir.size(); ++k) {
            std::string n = ofFilePath::getFileName(gdir.getPath(k));
            if (n.find("transcript") != std::string::npos || n.find("region") != std::string::npos) { txPath = gdir.getPath(k); break; }
        }
    }
    if (txPath.empty()) ofLogWarning("ofApp") << "no transcript-regions BED in " << dsDir << "genome/ (target track off for '" << d.name << "')";
    loadTxRegions(txPath);

    // BAM for GW phase/DC: datasets.json "bam", else a *.bam in <dataset>/bam/ or <dataset>/
    bamPath = d.bam;
    if (bamPath.empty()) {
        for (const std::string& sub : { std::string("bam"), std::string(".") }) {
            ofDirectory bdir(dsDir + sub); bdir.allowExt("bam"); bdir.listDir();
            if (bdir.size() > 0) { bamPath = bdir.getPath(0); break; }
        }
    }
    ofLogNotice("ofApp") << "GW BAM for '" << d.name << "': " << (bamPath.empty() ? "(none -> phase/DC disabled)" : bamPath);

    recomputeLayout();
}

// ---- coverage ----------------------------------------------------------------

void ofApp::openReaders() {
    readers.clear();
    for (const auto& p : bwPaths) {
        auto r = std::unique_ptr<BwReader>(new BwReader(p));
        if (!r->ok()) ofLogError("ofApp") << "cannot open bigWig: " << p;
        readers.push_back(std::move(r));
    }
}

std::vector<float> ofApp::sumBinnedChrom(const std::string& chrom, int nBins) {
    std::vector<float> acc(nBins, 0.f);
    for (auto& r : readers) { if (!r->ok()) continue;
        auto v = r->binnedChrom(chrom, nBins);
        for (int i = 0; i < nBins && i < (int)v.size(); ++i) acc[i] += v[i]; }
    return acc;
}

std::vector<float> ofApp::sumBinned(const std::string& chrom, long long s, long long e, int nBins) {
    std::vector<float> acc(nBins, 0.f);
    for (auto& r : readers) { if (!r->ok()) continue;
        auto v = r->binned(chrom, (uint32_t)s, (uint32_t)e, (uint32_t)nBins);
        for (int i = 0; i < nBins && i < (int)v.size(); ++i) acc[i] += v[i]; }
    return acc;
}

void ofApp::loadChroms() {
    arcs.clear();
    if (readers.empty() || !readers[0]->ok()) return;
    std::vector<BwReader::Chrom> ch;
    auto low = [](char c){ return (c >= 'A' && c <= 'Z') ? char(c + 32) : c; };
    // placed chromosomes: up to 6 leading letters then "ch"[r] + a digit — matches chr11_1, Chr01,
    // SPIMPch01; drops scaffolds/contigs (unitig.., ontctg.., scaffold_..) which lack it near the start.
    auto isChromLike = [&](const std::string& n) -> bool {
        for (std::size_t i = 0; i < n.size() && i <= 6 &&
             ((n[i] >= 'A' && n[i] <= 'Z') || (n[i] >= 'a' && n[i] <= 'z')); ++i) {
            if (low(n[i]) == 'c' && i + 1 < n.size() && low(n[i + 1]) == 'h') {
                std::size_t j = i + 2;
                if (j < n.size() && low(n[j]) == 'r') ++j;
                if (j < n.size() && n[j] >= '0' && n[j] <= '9') return true;
            }
        }
        return false;
    };
    for (auto& c : readers[0]->chroms())
        if (isChromLike(c.name)) ch.push_back(c);
    std::sort(ch.begin(), ch.end(),
              [](const BwReader::Chrom& a, const BwReader::Chrom& b) { return a.length > b.length; });
    int n = std::min((int)ch.size(), 24);
    for (int i = 0; i < n; ++i) { Arc a; a.id = ch[i].name; a.size = (long long)ch[i].length; arcs.push_back(a); }
}

void ofApp::layoutArcs() {
    long long total = 0; for (const auto& a : arcs) total += a.size;
    if (total == 0) return;
    const float gap = 0.012f;
    float avail = TWO_PI - gap * arcs.size();
    float ang = -HALF_PI;
    for (size_t i = 0; i < arcs.size(); ++i) {
        float span = avail * (float)((double)arcs[i].size / (double)total);
        arcs[i].a0 = ang; arcs[i].a1 = ang + span;
        arcs[i].color = ofColor::fromHsb((int)(i * 41) % 255, 170, 235);
        ang += span + gap;
    }
}

void ofApp::loadTrack(int binsPerArc) {
    trackBins.assign(arcs.size(), {}); trackMax = 1.f;
    for (size_t i = 0; i < arcs.size(); ++i) {
        auto v = sumBinnedChrom(arcs[i].id, binsPerArc);
        for (float x : v) trackMax = std::max(trackMax, x);
        trackBins[i] = std::move(v);
    }
}

long long ofApp::arcSize(const std::string& chrom) const {
    for (const auto& a : arcs) if (a.id == chrom) return a.size;
    return 0;
}

// ---- data --------------------------------------------------------------------

void ofApp::loadTsv(const std::string& path) {
    rows.clear(); rowsByChrom.clear(); rowsByTrigger.clear();
    std::ifstream f(path);
    if (!f) { ofLogError("ofApp") << "no tsv: " << path; return; }
    std::string line;
    if (!std::getline(f, line)) return;
    auto h = splitTab(line);
    int cChrom=-1,cSite=-1,cStrand=-1,cTarget=-1,cTrig=-1,cB=-1,cRS=-1,cRE=-1,cAllen=-1,cLeft=-1;
    for (int i = 0; i < (int)h.size(); ++i) {
        if (h[i]=="chrom") cChrom=i; else if (h[i]=="site") cSite=i;
        else if (h[i]=="strand") cStrand=i; else if (h[i]=="target") cTarget=i;
        else if (h[i]=="trigger") cTrig=i; else if (h[i]=="bscore") cB=i;
        else if (h[i]=="region_start") cRS=i; else if (h[i]=="region_end") cRE=i;
        else if (h[i]=="allen") cAllen=i; else if (h[i]=="left_of_sRNA") cLeft=i;
    }
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto v = splitTab(line);
        if ((int)v.size() <= cB) continue;
        Row r;
        r.chrom=v[cChrom]; r.site=atoll(v[cSite].c_str());
        r.strand = cStrand>=0 ? v[cStrand] : ".";
        r.target=v[cTarget]; r.trigger=v[cTrig]; r.bscore=(float)atof(v[cB].c_str());
        r.targetLower=lc(r.target); r.triggerLower=lc(r.trigger);
        if (cAllen>=0 && cAllen<(int)v.size()) r.allen=v[cAllen];
        if (cLeft>=0  && cLeft <(int)v.size()) r.leftOfSRNA=v[cLeft];
        if (cRS>=0 && cRS<(int)v.size()) r.regionStart=atoll(v[cRS].c_str());
        if (cRE>=0 && cRE<(int)v.size()) r.regionEnd=atoll(v[cRE].c_str());
        int idx = (int)rows.size();
        rowsByChrom[r.chrom].push_back(idx);
        rowsByTrigger[r.trigger].push_back(idx);
        rows.push_back(std::move(r));
    }
    ofLogNotice("ofApp") << "tsv rows: " << rows.size();
}

void ofApp::loadOrigins(const std::string& path) {
    origins.clear();
    std::ifstream f(path);
    if (!f) { ofLogWarning("ofApp") << "no origins: " << path; return; }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto v = splitTab(line);
        if (v.size() < 4) continue;
        Origin o; o.chrom=v[0]; o.pos=atoll(v[1].c_str());
        o.end    = (v.size() > 2) ? atoll(v[2].c_str()) : o.pos;
        o.strand = (v.size() > 4) ? v[4] : ".";
        origins[v[3]].push_back(o);
    }
    ofLogNotice("ofApp") << "origin triggers: " << origins.size();
}

void ofApp::loadTxRegions(const std::string& path) {
    txRegions.clear();
    if (path.empty()) return;
    std::ifstream f(path);
    if (!f) { ofLogWarning("ofApp") << "no transcript regions: " << path; return; }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        auto v = splitTab(line);
        if (v.size() < 3) continue;
        TxFeat t; t.start = atoll(v[1].c_str()); t.end = atoll(v[2].c_str());
        t.name   = (v.size() > 3) ? v[3] : "";
        t.strand = (v.size() > 5) ? v[5] : ((v.size() > 4) ? v[4] : ".");
        txRegions[v[0]].push_back(t);
    }
    size_t n = 0; for (auto& kv : txRegions) n += kv.second.size();
    ofLogNotice("ofApp") << "transcript regions: " << n << " on " << txRegions.size() << " chroms";
}

void ofApp::loadFasta(const std::string& path) {
    triggerSeq.clear();
    std::ifstream f(path);
    if (!f) { ofLogWarning("ofApp") << "no fasta: " << path; return; }
    std::string line, id;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();   // tolerate CRLF
        if (line.empty()) continue;
        if (line[0] == '>') {
            id = line.substr(1);
            // Key by the FULL header (trailing whitespace trimmed) so miR names that
            // contain spaces ("nta-miR6162 MIMAT0024777 Nicotiana tabacum miR6162")
            // match the tsv 'trigger' field verbatim. sRNA names have no spaces, so
            // the full header == the old first-token key -> they are unaffected.
            while (!id.empty() && (id.back() == ' ' || id.back() == '\t')) id.pop_back();
            triggerSeq[id] = "";
        } else if (!id.empty()) {
            triggerSeq[id] += line;                  // sequence (append multi-line)
        }
    }
    ofLogNotice("ofApp") << "fasta seqs: " << triggerSeq.size();
}

void ofApp::loadClusterTx(const std::string& path) {
    clusterTx.clear();
    if (path.empty()) return;
    std::ifstream f(path);
    if (!f) { ofLogWarning("ofApp") << "no clusters_on_tx: " << path; return; }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto v = splitTab(line);
        if (v.size() < 4) continue;
        // name = transcript|Cluster_N|DC=..|MIR=..|PT=..|AX=..
        std::string name = v[3];
        size_t p1 = name.find('|');
        if (p1 == std::string::npos) continue;
        std::string tx = name.substr(0, p1);
        size_t p2 = name.find('|', p1 + 1);
        std::string cl = name.substr(p1 + 1, (p2 == std::string::npos ? name.size() : p2) - p1 - 1);
        std::string& cur = clusterTx[cl];
        if (cur.empty()) cur = tx;                       // join isoforms if a cluster maps to several
        else if (cur.find(tx) == std::string::npos) cur += "," + tx;
    }
    clusterTxLower.clear();
    for (auto& kv : clusterTx) clusterTxLower[kv.first] = lc(kv.second);
    ofLogNotice("ofApp") << "cluster->tx: " << clusterTx.size();
}

void ofApp::filterTable(const std::string& chrom) {
    activeChrom=chrom; tableMode = TableMode::ByChrom; tableRows.clear();
    auto it = rowsByChrom.find(chrom);
    if (it != rowsByChrom.end()) tableRows = it->second;
    std::sort(tableRows.begin(), tableRows.end(),
              [this](int a, int b){ return rows[a].bscore > rows[b].bscore; });
    tableScroll=0; selRow=-1; selOrigins.clear();
    applySearch();
}

void ofApp::showTargetsOfTrigger(const std::string& trig) {
    tableMode = TableMode::ByTrigger; triggerFilter = trig; tableRows.clear();
    auto it = rowsByTrigger.find(trig);
    if (it != rowsByTrigger.end()) tableRows = it->second;
    std::sort(tableRows.begin(), tableRows.end(),
              [this](int a, int b){ return rows[a].bscore > rows[b].bscore; });
    tableScroll = 0;                     // keep selRow selected
    applySearch();
}

void ofApp::backToChromTable() {
    tableMode = TableMode::ByChrom; tableRows.clear();
    auto it = rowsByChrom.find(activeChrom);
    if (it != rowsByChrom.end()) tableRows = it->second;
    std::sort(tableRows.begin(), tableRows.end(),
              [this](int a, int b){ return rows[a].bscore > rows[b].bscore; });
    tableScroll = 0;                     // keep selRow selected
    applySearch();
}

void ofApp::applySearch() {
    viewRows.clear();
    if (searchQuery.empty()) {
        viewRows = tableRows;
    } else {
        std::string q = lc(searchQuery);              // lowercase the query ONCE
        for (int idx : tableRows) {
            if (rows[idx].targetLower.find(q)  != std::string::npos) { viewRows.push_back(idx); continue; }
            if (rows[idx].triggerLower.find(q) != std::string::npos) { viewRows.push_back(idx); continue; }  // search trigger too
            auto it = clusterTxLower.find(rows[idx].target);   // precomputed lower transcript
            if (it != clusterTxLower.end() && it->second.find(q) != std::string::npos)
                viewRows.push_back(idx);
        }
    }
    tableScroll = 0; copyLo = copyHi = -1;
}

void ofApp::copySelection() {
    if (copyLo < 0 || viewRows.empty()) return;
    int lo = std::max(0, std::min(copyLo, copyHi));
    int hi = std::min((int)viewRows.size() - 1, std::max(copyLo, copyHi));
    std::string out;
    for (int i = lo; i <= hi; ++i) {
        const Row& r = rows[viewRows[i]];
        std::string tx, seq;
        auto it = clusterTx.find(r.target);   if (it != clusterTx.end()) tx = it->second;
        auto st = triggerSeq.find(r.trigger); if (st != triggerSeq.end()) seq = st->second;
        out += r.trigger + "\t" + r.target + "\t" + tx + "\t" + r.chrom + "\t"
             + ofToString(r.site) + "\t" + r.allen + "\t" + r.leftOfSRNA + "\t"
             + ofToString(r.bscore, 4) + "\t" + seq + "\n";
    }
    ofSetClipboardString(out);
    toast = "copied " + ofToString(hi - lo + 1) + " row(s)";
    toastUntil = ofGetElapsedTimef() + 1.5f;
}

void ofApp::selectRow(int idx) {
    selRow = idx; selOrigins.clear();
    if (idx < 0) return;
    auto it = origins.find(rows[idx].trigger);
    if (it == origins.end()) return;
    for (const auto& o : it->second)
        if (arcSize(o.chrom) > 0) selOrigins.push_back(o);   // keep origins on drawn arcs
}

// ---- focus (zoom on one origin) ----------------------------------------------

void ofApp::enterFocus(int oi) {
    if (oi < 0 || oi >= (int)selOrigins.size()) return;
    focusOrigin = selOrigins[oi];
    focusTarget = false;
    long long L = arcSize(focusOrigin.chrom);
    viewStart = std::max(0LL, focusOrigin.pos - originPad);
    viewEnd = focusOrigin.pos + originPad;
    if (L > 0) viewEnd = std::min(viewEnd, L);
    mode = Mode::Focus;
    rebinFocus();
}

void ofApp::enterFocusTarget() {
    if (selRow < 0) return;
    const Row& r = rows[selRow];
    focusOrigin = { r.chrom, r.site };           // reuse focus machinery, pointed at the target
    focusTarget = true;
    long long L = arcSize(r.chrom);
    long long s = r.regionStart, e = r.regionEnd;
    if (e <= s) { s = std::max(0LL, r.site - originPad); e = r.site + originPad; }
    if (L > 0) e = std::min(e, L);
    viewStart = std::max(0LL, s); viewEnd = e;
    mode = Mode::Focus;
    rebinFocus();
}

bool ofApp::targetDotHit(int x, int y) const {
    if (selRow < 0) return false;
    ofPoint p;
    if (!genomicToPoint(rows[selRow].chrom, rows[selRow].site, p)) return false;
    return ofDist(x, y, center.x + p.x, center.y + p.y) < 12.f;
}

void ofApp::rebinFocus() {
    focusMax = 1.f;
    focusBins = sumBinned(focusOrigin.chrom, viewStart, viewEnd, focusResBins);
    for (float v : focusBins) focusMax = std::max(focusMax, v);
}

// ---- geometry ----------------------------------------------------------------

void ofApp::recomputeLayout() {
    center = ofPoint(ofGetWidth() * 0.5f, layoutTop() * 0.5f);
    radius = std::min((float)ofGetWidth(), layoutTop()) * 0.36f;
}

float ofApp::angleAt(int x, int y) const { return atan2f((float)y - center.y, (float)x - center.x); }

int ofApp::arcAtAngle(float ang) const {
    auto norm = [](float a){ a = fmodf(a, TWO_PI); return a < 0 ? a + TWO_PI : a; };
    float q = norm(ang);
    for (size_t i = 0; i < arcs.size(); ++i) {
        float a0 = norm(arcs[i].a0 + rotation), a1 = norm(arcs[i].a1 + rotation);
        bool inside = (a0 <= a1) ? (q >= a0 && q <= a1) : (q >= a0 || q <= a1);
        if (inside) return (int)i;
    }
    return -1;
}

bool ofApp::genomicToPoint(const std::string& chrom, long long pos, ofPoint& out) const {
    for (const auto& a : arcs) {
        if (a.id != chrom) continue;
        float frac = a.size > 0 ? ofClamp((float)((double)pos/(double)a.size), 0.f, 1.f) : 0.f;
        float t = ofLerp(a.a0 + rotation, a.a1 + rotation, frac);
        out.set((radius - 4)*cosf(t), (radius - 4)*sinf(t));
        return true;
    }
    return false;
}

bool ofApp::originAngle(const Origin& o, float& ang) const {
    for (const auto& a : arcs) {
        if (a.id != o.chrom) continue;
        float frac = a.size > 0 ? ofClamp((float)((double)o.pos/(double)a.size), 0.f, 1.f) : 0.f;
        ang = ofLerp(a.a0 + rotation, a.a1 + rotation, frac);
        return true;
    }
    return false;
}

bool ofApp::focusArrowAngle(float& ang) const {
    long long span = viewEnd - viewStart;
    if (span <= 0) return false;
    float frac = ofClamp((float)((double)(focusOrigin.pos - viewStart)/(double)span), 0.f, 1.f);
    ang = ofLerp(fa0, fa1, frac);
    return true;
}

int ofApp::arrowAtPixel(int x, int y) const {
    for (size_t i = 0; i < selOrigins.size(); ++i) {
        float ang;
        if (!originAngle(selOrigins[i], ang)) continue;
        ofPoint p(center.x + (radius + arrowR)*cosf(ang), center.y + (radius + arrowR)*sinf(ang));
        if (ofDist(x, y, p.x, p.y) < 16.f) return (int)i;
    }
    return -1;
}

// ---- draw --------------------------------------------------------------------

void ofApp::draw() {
    ofBackground(18);
    ofPushMatrix();
    ofTranslate(center);
    if (mode == Mode::Overview) { drawRing(); drawChords(); }
    else                          drawFocus();
    drawArrows();
    ofPopMatrix();
    drawTable();
    gui.draw();

    ofSetColor(230);
    if (mode == Mode::Focus) {
        long long oe = (focusOrigin.end > focusOrigin.pos) ? focusOrigin.end : focusOrigin.pos;
        std::string coord = focusOrigin.chrom + ":" + ofToString(focusOrigin.pos)
                          + (oe > focusOrigin.pos ? "-" + ofToString(oe) : "")
                          + (focusOrigin.strand.empty() || focusOrigin.strand == "." ? "" : "(" + focusOrigin.strand + ")");
        std::string head = std::string(focusTarget ? "target  " : "trigger origin  ") + coord
                           + "   [" + ofToString(viewStart) + "-" + ofToString(viewEnd) + "]";
        ofDrawBitmapString(head, 12, 18);
        float bx = 12 + head.size() * 8.f + 8.f;            // just after the ']'
        gwBtn.set(bx, 6, 28, 15);
        ofSetColor(58, 110, 88); ofDrawRectangle(gwBtn);
        ofSetColor(225, 255, 235); ofDrawBitmapString("GW", bx + 6, 18);
        ofSetColor(160);
        ofDrawBitmapString("scroll: zoom   Esc: back", bx + 36, 18);
    }
    else if (selRow >= 0)
        ofDrawBitmapString(ofToString((int)selOrigins.size())
                           + " origin(s) — red arrow: zoom trigger origin   yellow dot: zoom target", 12, 18);

    // trigger sequence for the selected row
    if (selRow >= 0) {
        auto it = triggerSeq.find(rows[selRow].trigger);
        std::string seq = (it != triggerSeq.end()) ? it->second : "(not in fasta)";
        ofSetColor(120, 255, 160);
        ofDrawBitmapString(rows[selRow].trigger + "  " + seq + "  (" + ofToString((int)seq.size()) + " nt)",
                           12, layoutTop() - 10);
    }

    if (ctxOpen) {                                   // right-click context menu (on top)
        ofSetColor(64); ofDrawRectangle(ctxX, ctxY, 300, 20);
        ofSetColor(235);
        ofDrawBitmapString("show other targets of this trigger", ctxX + 6, ctxY + 14);
    }

    if (gwOpen) drawGW();                            // GW overlay sits above everything
}

void ofApp::drawRing() {
    for (size_t i = 0; i < arcs.size(); ++i) {
        const Arc& a = arcs[i];
        float a0 = a.a0 + rotation, a1 = a.a1 + rotation;
        bool sel = ((int)i == selectedArc);
        ofSetColor(sel ? ofColor::yellow : a.color);
        ofSetLineWidth(sel ? 3.f : 1.5f);
        ofPolyline pl;
        int steps = std::max(2, (int)((a1 - a0) / 0.02f));
        for (int s = 0; s <= steps; ++s) { float t = ofLerp(a0, a1, (float)s/steps);
            pl.addVertex(radius*cosf(t), radius*sinf(t)); }
        pl.draw();
        if (i < trackBins.size() && !trackBins[i].empty()) {
            const auto& bins = trackBins[i];
            ofSetColor(a.color, 190); ofSetLineWidth(1.f);
            for (size_t b = 0; b < bins.size(); ++b) {
                float t = ofLerp(a0, a1, (b + 0.5f)/bins.size());
                float h = (bins[b]/trackMax)*trackDepth;
                ofDrawLine((radius-2)*cosf(t),(radius-2)*sinf(t),(radius-2-h)*cosf(t),(radius-2-h)*sinf(t));
            }
        }
        float tm = (a0 + a1)*0.5f;
        ofSetColor(sel ? ofColor::yellow : ofColor(200));
        ofDrawBitmapString(a.id, (radius+8)*cosf(tm), (radius+8)*sinf(tm));
    }
}

void ofApp::drawChords() {
    if (selRow < 0) return;
    const Row& r = rows[selRow];
    ofPoint sitePt;
    bool haveSite = genomicToPoint(r.chrom, r.site, sitePt);
    ofSetColor(120, 200, 255, 160); ofSetLineWidth(1.5f);
    for (const auto& o : selOrigins) {
        ofPoint op;
        if (!genomicToPoint(o.chrom, o.pos, op) || !haveSite) continue;
        ofPolyline chord;
        for (int s = 0; s <= 40; ++s) { float t = s/40.f, u = 1-t;
            chord.addVertex(u*u*op.x + t*t*sitePt.x, u*u*op.y + t*t*sitePt.y); }
        chord.draw();
    }
    if (haveSite) { ofSetColor(255, 220, 60); ofDrawCircle(sitePt, 4); }  // target slice-site
}

void ofApp::drawArrows() {
    if (mode == Mode::Focus) {
        ofSetColor(focusTarget ? ofColor(255, 220, 60) : ofColor(255, 60, 60));
        float ang; if (focusArrowAngle(ang)) drawArrowAt(ang, radius, arrowR);
        return;
    }
    ofSetColor(255, 60, 60);                 // origin arrows
    for (const auto& o : selOrigins) {
        float ang; if (originAngle(o, ang)) drawArrowAt(ang, radius, arrowR);
    }
}

void ofApp::drawFocus() {
    ofSetColor(180); ofSetLineWidth(3.f);
    ofPolyline pl;
    for (int s = 0; s <= 240; ++s) { float t = ofLerp(fa0, fa1, (float)s/240);
        pl.addVertex(radius*cosf(t), radius*sinf(t)); }
    pl.draw();
    if (!focusBins.empty()) {
        ofSetColor(120, 200, 255, 210); ofSetLineWidth(1.f);
        for (size_t b = 0; b < focusBins.size(); ++b) {
            float t = ofLerp(fa0, fa1, (b + 0.5f)/focusBins.size());
            float h = (focusBins[b]/focusMax)*trackDepth*1.8f;
            ofDrawLine((radius-2)*cosf(t),(radius-2)*sinf(t),(radius-2-h)*cosf(t),(radius-2-h)*sinf(t));
        }
    }
    ofSetColor(200);
    ofDrawBitmapString(ofToString(viewStart), (radius+8)*cosf(fa0), (radius+8)*sinf(fa0));
    ofDrawBitmapString(ofToString(viewEnd),   (radius+8)*cosf(fa1), (radius+8)*sinf(fa1));
}

void ofApp::drawTable() {
    float top = layoutTop();
    ofSetColor(28); ofDrawRectangle(0, top, ofGetWidth(), ofGetHeight() - top);

    // search bar
    ofSetColor(searchFocused ? ofColor(45, 62, 82) : ofColor(42));
    ofDrawRectangle(12, top + 2, ofGetWidth() - 24, 18);
    ofSetColor(searchFocused ? 255 : 175);
    ofDrawBitmapString("search target/trigger: " + searchQuery + (searchFocused ? "_" : ""), 18, top + 15);

    // header
    ofSetColor(160);
    std::string hdr;
    if (tableMode == TableMode::ByTrigger)
        hdr = "targets of " + triggerFilter + "   —   " + ofToString((int)viewRows.size()) + " rows";
    else
        hdr = activeChrom.empty()
            ? "click a chromosome arc to list its triggers"
            : (activeChrom + "  —  " + ofToString((int)viewRows.size())
               + " triggers   (drag: select rows   c: copy   right-click trigger: other targets)");
    ofDrawBitmapString(hdr, 12, top + 38);
    if (tableMode == TableMode::ByTrigger) {          // back button
        ofSetColor(90); ofDrawRectangle(ofGetWidth() - 96, top + 24, 84, 18);
        ofSetColor(235); ofDrawBitmapString("< back", ofGetWidth() - 88, top + 37);
    }

    // rows (from the search-filtered viewRows)
    float y0 = top + 56;
    int visible = (int)((ofGetHeight() - y0) / rowH);
    for (int i = tableScroll; i < tableScroll + visible && i < (int)viewRows.size(); ++i) {
        int idx = viewRows[i]; const Row& r = rows[idx];
        float yy = y0 + (i - tableScroll) * rowH;
        if (idx == selRow) { ofSetColor(60, 90, 120); ofDrawRectangle(0, yy-11, ofGetWidth(), rowH); }
        ofSetColor(idx == selRow ? ofColor(255) : ofColor(210));
        std::string tgt = r.target;
        auto tit = clusterTx.find(r.target);
        if (tit != clusterTx.end()) tgt += " [" + tit->second + "]";   // Soltu.DM transcript
        ofDrawBitmapString(r.trigger + "   " + tgt + "   " + ofToString(r.site)
                           + "   allen=" + (r.allen.empty()?"NA":r.allen)
                           + "   L=" + (r.leftOfSRNA.empty()?"NA":r.leftOfSRNA)
                           + "   b=" + ofToString(r.bscore, 4), 12, yy);
    }

    if (ofGetElapsedTimef() < toastUntil) {           // "copied" feedback
        ofSetColor(120, 255, 160);
        ofDrawBitmapString(toast, 12, ofGetHeight() - 6);
    }
}

// ---- interaction -------------------------------------------------------------

// ---- GW genome-window view ---------------------------------------------------

void ofApp::openGW() {
    if (selRow < 0) return;
    // bottom pane = the current zoom window, exactly as shown in [start-end]
    gwBotChrom = focusOrigin.chrom; gwBotS = viewStart; gwBotE = viewEnd;
    gwBotMark  = focusOrigin.pos;                       // focal locus (slice site / origin)
    // top pane = the locus the trigger comes from (its genomic origin), +/- originPad of context
    Origin org = (focusTarget && !selOrigins.empty()) ? selOrigins[0] : focusOrigin;
    long long L = arcSize(org.chrom);
    gwTopChrom = org.chrom;
    gwTopS = std::max(0LL, org.pos - originPad);
    gwTopE = org.pos + originPad;
    if (L > 0) gwTopE = std::min(gwTopE, L);
    gwTopMark = org.pos;
    gwTopMarkEnd = (org.end > org.pos) ? org.end : org.pos;
    gwTopStrand = (org.strand.empty() || org.strand == ".") ? "" : "(" + org.strand + ")";
    // coverage for both windows (fixed resolution; mapped to pane width when drawn)
    gwTopBins = sumBinned(gwTopChrom, gwTopS, gwTopE, gwBinN);
    gwBotBins = sumBinned(gwBotChrom, gwBotS, gwBotE, gwBinN);
    gwTopMax = 1.f; for (float v : gwTopBins) gwTopMax = std::max(gwTopMax, v);
    gwBotMax = 1.f; for (float v : gwBotBins) gwBotMax = std::max(gwBotMax, v);
    gwTrigger = rows[selRow].trigger;
    gwRevealed.clear();
    gwOpen = true; gwDragSplit = false;
    requestGwPhase();
}

int ofApp::gwSliceHit(const ofRectangle& pane, const std::string& chrom,
                      long long s, long long e, int x, int y) const {
    if (e <= s) return -1;
    float plotL = pane.x + 8, plotR = pane.getRight() - 8, plotW = plotR - plotL;
    float topY = pane.y + 18, baseY = pane.getBottom() - 22;
    if (y < topY || y > baseY || x < plotL - 3 || x > plotR + 3) return -1;
    auto it = rowsByChrom.find(chrom);
    if (it == rowsByChrom.end()) return -1;
    int best = -1; float bestd = 4.f;                          // 4 px tolerance
    for (int ri : it->second) {
        long long st = rows[ri].site;
        if (st < s || st > e) continue;
        float sx = plotL + (float)((double)(st - s) / (double)(e - s)) * plotW;
        float d = std::fabs(sx - (float)x);
        if (d < bestd) { bestd = d; best = ri; }
    }
    return best;
}

void ofApp::requestGwPhase() {
    if (bamPath.empty() || gwScript.empty()) {
        std::lock_guard<std::mutex> lk(gwPhaseMtx); gwTopPhase = ""; gwBotPhase = ""; return;
    }
    { std::lock_guard<std::mutex> lk(gwPhaseMtx); gwTopPhase = "phase: computing..."; gwBotPhase = "phase: computing..."; }
    gwPhasePending = true;
    gwPhaseChangeT = ofGetElapsedTimef();
}

std::string ofApp::queryGwPhase(const std::string& chrom, long long s, long long e) {
    if (bamPath.empty() || gwScript.empty()) return "";
    std::string cmd = gwPython + " '" + gwScript + "' --bam '" + bamPath + "' --chrom '" + chrom
                    + "' --start " + std::to_string(s) + " --end " + std::to_string(e) + " 2>/dev/null";
    FILE* f = popen(cmd.c_str(), "r");
    if (!f) return "";
    char buf[256]; std::string out;
    if (fgets(buf, sizeof(buf), f)) out = buf;
    pclose(f);
    while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) out.pop_back();
    return out;
}

void ofApp::update() {
    // debounce: recompute GW phase/DC once the window has been stable ~0.35 s
    if (gwOpen && gwPhasePending && !gwPhaseBusy.load()
        && (ofGetElapsedTimef() - gwPhaseChangeT) > 0.35f) {
        gwPhasePending = false;
        gwPhaseBusy = true;
        std::string tc = gwTopChrom, bc = gwBotChrom;         // snapshot for the worker
        long long ts = gwTopS, te = gwTopE, bs = gwBotS, be = gwBotE;
        std::thread([this, tc, ts, te, bc, bs, be]() {
            std::string t = queryGwPhase(tc, ts, te);
            std::string b = queryGwPhase(bc, bs, be);
            { std::lock_guard<std::mutex> lk(gwPhaseMtx); gwTopPhase = t; gwBotPhase = b; }
            gwPhaseBusy = false;
        }).detach();
    }
}

void ofApp::openJBrowse(const std::string& chrom, long long s, long long e) {
    if (jbUrl.empty()) return;
    std::string sep = (jbUrl.find('?') != std::string::npos) ? "&" : "?";
    ofLaunchBrowser(jbUrl + sep + "loc=" + chrom + ":" + ofToString(s) + "-" + ofToString(e));
}

void ofApp::drawCovPane(ofRectangle rc, const std::string& title, const std::string& chrom,
                        long long s, long long e, long long mark,
                        const std::vector<float>& bins, float mx, ofColor col) {
    ofSetColor(15); ofDrawRectangle(rc);
    ofSetColor(60); ofNoFill(); ofDrawRectangle(rc); ofFill();
    const float padL = 8, padR = 8, padTop = 18, padBot = 22;   // padBot leaves a BED feature lane
    float plotL = rc.x + padL, plotR = rc.getRight() - padR;
    float baseY = rc.getBottom() - padBot, topY = rc.y + padTop;
    float plotW = plotR - plotL, plotH = std::max(4.f, baseY - topY);

    ofSetColor(210);
    ofDrawBitmapString(title, rc.x + 6, rc.y + 13);

    int N = (int)bins.size();
    if (N > 0 && mx > 0.f) {
        ofSetColor(col); ofSetLineWidth(1.f);
        for (int b = 0; b < N; ++b) {
            float px = plotL + (b + 0.5f) / N * plotW;
            float h  = bins[b] / mx * plotH;
            ofDrawLine(px, baseY, px, baseY - h);
        }
    }
    ofSetColor(90); ofDrawLine(plotL, baseY, plotR, baseY);          // baseline

    // predicted slice sites in this window (from the bscore TSV / --slice_bed): all darker yellow...
    auto rit = rowsByChrom.find(chrom);
    if (rit != rowsByChrom.end() && e > s) {
        ofSetLineWidth(1.f); ofSetColor(150, 130, 30);
        for (int ri : rit->second) {
            long long st = rows[ri].site;
            if (st < s || st > e || st == mark) continue;
            float sx = plotL + (float)((double)(st - s) / (double)(e - s)) * plotW;
            ofDrawLine(sx, topY, sx, baseY);
        }
    }
    // ...and the selected/focal slice site in bright yellow, on top
    if (e > s && mark >= s && mark <= e) {
        float mxp = plotL + (float)((double)(mark - s) / (double)(e - s)) * plotW;
        ofSetColor(255, 220, 60); ofSetLineWidth(1.5f);
        ofDrawLine(mxp, topY, mxp, baseY);
        ofSetLineWidth(1.f);
        ofDrawBitmapString(ofToString(mark), mxp + 3, topY + 10);
    }

    // labels for revealed dim slice sites (click a line to toggle): trigger + its origin
    for (int ri : gwRevealed) {
        if (ri < 0 || ri >= (int)rows.size() || rows[ri].chrom != chrom) continue;
        long long st = rows[ri].site;
        if (st < s || st > e) continue;
        float sx = plotL + (float)((double)(st - s) / (double)(e - s)) * plotW;
        std::string lbl = rows[ri].trigger;
        auto oi = origins.find(rows[ri].trigger);
        if (oi != origins.end() && !oi->second.empty())
            lbl += "  @" + oi->second[0].chrom + ":" + ofToString(oi->second[0].pos);
        ofSetColor(255, 235, 130); ofSetLineWidth(1.5f);
        ofDrawLine(sx, topY, sx, baseY);                    // re-draw this line brighter while labeled
        ofSetLineWidth(1.f);
        ofDrawBitmapString(lbl, sx + 3, topY + 22 + (ri % 6) * 11);   // stagger to reduce overlap
    }

    // transcript-region features overlapping the window (BED track lane just below the baseline)
    auto it = txRegions.find(chrom);
    if (it != txRegions.end() && e > s) {
        float laneY = baseY + 3, laneH = 7.f;
        for (const auto& t : it->second) {
            if (t.end < s || t.start > e) continue;
            float x0 = plotL + (float)((double)(std::max(t.start, s) - s) / (double)(e - s)) * plotW;
            float x1 = plotL + (float)((double)(std::min(t.end,   e) - s) / (double)(e - s)) * plotW;
            if (x1 - x0 < 1.5f) x1 = x0 + 1.5f;
            ofSetColor(90, 200, 120); ofDrawRectangle(x0, laneY, x1 - x0, laneH);
            if (x1 - x0 > 52.f && !t.name.empty()) {                 // label only if wide enough
                ofSetColor(210, 255, 220);
                ofDrawBitmapString(t.name, x0 + 2, laneY + laneH + 8);
            }
        }
    }

    ofSetColor(140);
    ofDrawBitmapString(ofToString(s), plotL, rc.getBottom() - 3);
    std::string es = ofToString(e);
    ofDrawBitmapString(es, plotR - es.size() * 8.f, rc.getBottom() - 3);
}

void ofApp::drawGW() {
    float W = ofGetWidth(), H = ofGetHeight();
    ofSetColor(0, 0, 0, 190); ofDrawRectangle(0, 0, W, H);          // dim background
    ofRectangle panel(50, 40, W - 100, H - 80);
    ofSetColor(24); ofDrawRectangle(panel);
    ofSetColor(80); ofNoFill(); ofDrawRectangle(panel); ofFill();

    ofSetColor(225);
    ofDrawBitmapString("Genome window   trigger " + (selRow >= 0 ? rows[selRow].trigger : std::string()),
                       panel.x + 10, panel.y + 16);
    gwCloseBtn.set(panel.getRight() - 26, panel.y + 4, 18, 16);      // close
    ofSetColor(120, 50, 50); ofDrawRectangle(gwCloseBtn);
    ofSetColor(255); ofDrawBitmapString("x", gwCloseBtn.x + 6, gwCloseBtn.y + 13);

    float cTop = panel.y + 26, cH = panel.getBottom() - 8 - cTop;
    gwContentTop = cTop; gwContentH = cH;
    const float splitH = 10.f;
    float topH = (cH - splitH) * ofClamp(gwSplit, 0.12f, 0.88f);
    float botH = (cH - splitH) - topH;
    gwTopPane.set(panel.x + 8, cTop, panel.width - 16, topH);
    gwBotPane.set(panel.x + 8, cTop + topH + splitH, panel.width - 16, botH);

    std::string topPhase, botPhase;
    { std::lock_guard<std::mutex> lk(gwPhaseMtx); topPhase = gwTopPhase; botPhase = gwBotPhase; }
    // top title = the trigger + its exact origin coordinate (from the BED), then phase/DC
    std::string topCoord = gwTopChrom + ":" + ofToString(gwTopMark)
                         + (gwTopMarkEnd > gwTopMark ? "-" + ofToString(gwTopMarkEnd) : "") + gwTopStrand;
    std::string topTitle = "trigger origin  " + gwTrigger + "  " + topCoord + "   " + topPhase;
    std::string botTitle = "window  " + gwBotChrom + ":" + ofToString(gwBotS) + "-" + ofToString(gwBotE)
                         + "   " + botPhase;
    drawCovPane(gwTopPane, topTitle, gwTopChrom, gwTopS, gwTopE, gwTopMark, gwTopBins, gwTopMax, ofColor(255, 120, 120));
    drawCovPane(gwBotPane, botTitle, gwBotChrom, gwBotS, gwBotE, gwBotMark, gwBotBins, gwBotMax, ofColor(120, 200, 255));

    gwTopJB.set(0, 0, 0, 0); gwBotJB.set(0, 0, 0, 0);      // 'JB' -> open this window in JBrowse
    if (!jbUrl.empty()) {
        auto jb = [&](ofRectangle& btn, const ofRectangle& pane) {
            btn.set(pane.getRight() - 26, pane.y + 3, 22, 14);
            ofSetColor(70, 90, 150); ofDrawRectangle(btn);
            ofSetColor(210, 225, 255); ofDrawBitmapString("JB", btn.x + 5, btn.y + 11);
        };
        jb(gwTopJB, gwTopPane); jb(gwBotJB, gwBotPane);
    }

    gwSplitBar.set(panel.x + 8, cTop + topH, panel.width - 16, splitH);
    ofSetColor(gwDragSplit ? ofColor(120, 180, 140) : ofColor(70)); ofDrawRectangle(gwSplitBar);
    ofSetColor(150);
    ofDrawBitmapString("= drag to resize =", gwSplitBar.getCenter().x - 72, gwSplitBar.y + 8);
}

void ofApp::mousePressed(int x, int y, int button) {
    // GW overlay is modal: it swallows all clicks while open
    if (gwOpen) {
        if (gwCloseBtn.inside(x, y)) { closeGW(); return; }
        if (gwTopJB.inside(x, y)) { openJBrowse(gwTopChrom, gwTopS, gwTopE); return; }
        if (gwBotJB.inside(x, y)) { openJBrowse(gwBotChrom, gwBotS, gwBotE); return; }
        if (gwSplitBar.inside(x, y)) { gwDragSplit = true; return; }
        int hit = gwSliceHit(gwTopPane, gwTopChrom, gwTopS, gwTopE, x, y);   // click a slice line -> toggle label
        if (hit < 0) hit = gwSliceHit(gwBotPane, gwBotChrom, gwBotS, gwBotE, x, y);
        if (hit >= 0) {
            auto p = std::find(gwRevealed.begin(), gwRevealed.end(), hit);
            if (p != gwRevealed.end()) gwRevealed.erase(p); else gwRevealed.push_back(hit);
        }
        return;
    }

    // let ofxGui/ofxDropdown own its area (panel + room for the expanded list)
    ofRectangle g = gui.getShape();
    if (ofRectangle(g.x - 4, g.y - 4, g.width + 8, g.height + datasets.size() * 26 + 44).inside(x, y))
        return;

    // "GW" button in the focus header -> open the genome-window view
    if (mode == Mode::Focus && button == OF_MOUSE_BUTTON_LEFT && gwBtn.inside(x, y)) { openGW(); return; }

    // right-click context menu: click the item to apply, anything else dismisses
    if (ctxOpen) {
        if (ofRectangle(ctxX, ctxY, 300, 20).inside(x, y)) showTargetsOfTrigger(ctxTrigger);
        ctxOpen = false; return;
    }
    float top = layoutTop();

    // search box focus toggle
    if (ofRectangle(12, top + 2, ofGetWidth() - 24, 18).inside(x, y)) { searchFocused = true; return; }
    searchFocused = false;

    // < back button (trigger-table mode)
    if (tableMode == TableMode::ByTrigger &&
        ofRectangle(ofGetWidth() - 96, top + 24, 84, 18).inside(x, y)) { backToChromTable(); return; }

    // right-click a row's first column (the trigger) -> context menu
    if (button == OF_MOUSE_BUTTON_RIGHT && y >= top + 56 && x < 280) {
        int i = tableScroll + (int)((y - (top + 56) + 11) / rowH);
        if (i >= 0 && i < (int)viewRows.size()) { ctxTrigger = rows[viewRows[i]].trigger; ctxOpen = true; ctxX = x; ctxY = y; }
        return;
    }
    if (button != OF_MOUSE_BUTTON_LEFT) return;      // the rest is left-click only

    // click the green trigger-sequence line (just above the table) -> copy as FASTA: >name \n seq
    if (selRow >= 0 && y >= top - 18 && y <= top - 2 && x >= 12) {
        const std::string& trig = rows[selRow].trigger;
        auto st = triggerSeq.find(trig);
        if (st != triggerSeq.end()) {
            ofSetClipboardString(">" + trig + "\n" + st->second);   // multi-origin: no single position
            toast = "copied FASTA"; toastUntil = ofGetElapsedTimef() + 1.5f;
        }
        return;
    }

    if (mode == Mode::Overview) {
        int ai = arrowAtPixel(x, y);
        if (ai >= 0) { enterFocus(ai); return; }         // red arrow -> zoom into origin
        if (targetDotHit(x, y)) { enterFocusTarget(); return; }  // yellow dot -> zoom into target
        if (y < top) {                                   // ring
            dragging = true; pressX = x; pressY = y;
            dragStartAngle = angleAt(x, y); rotAtDragStart = rotation;
            selectedArc = arcAtAngle(dragStartAngle);
        } else {                                         // table row: select
            int i = tableScroll + (int)((y - (top + 56) + 11) / rowH);
            if (i >= 0 && i < (int)viewRows.size()) selectRow(viewRows[i]);
        }
    } else {                                             // Focus: drag pans
        dragging = true; pressX = x; pressY = y;
        dragStartAngle = angleAt(x, y); viewStartAtDrag = viewStart;
    }
}

void ofApp::mouseDragged(int x, int y, int button) {
    if (gwOpen) {
        if (gwDragSplit && gwContentH > 1.f)
            gwSplit = ofClamp((y - gwContentTop) / gwContentH, 0.12f, 0.88f);
        return;
    }
    if (!dragging) return;
    if (mode == Mode::Focus) {
        long long span = viewEnd - viewStart, L = arcSize(focusOrigin.chrom);
        float dFrac = (angleAt(x, y) - dragStartAngle)/(fa1 - fa0);
        long long ns = viewStartAtDrag - (long long)(dFrac*span);
        ns = std::max(0LL, std::min(ns, L - span));
        viewStart = ns; viewEnd = ns + span; rebinFocus();
    } else {
        rotation = rotAtDragStart + (angleAt(x, y) - dragStartAngle);
    }
}

void ofApp::mouseReleased(int x, int y, int button) {
    if (gwOpen) { gwDragSplit = false; return; }
    if (mode == Mode::Overview && dragging && pressY < layoutTop() && ofDist(x, y, pressX, pressY) < 5.f) {
        int hit = arcAtAngle(angleAt(x, y));
        if (hit >= 0) { selectedArc = hit; filterTable(arcs[hit].id); }
    }
    dragging = false;
}

void ofApp::mouseScrolled(int x, int y, float sx, float sy) {
    if (gwOpen) {
        auto zoomPane = [&](const ofRectangle& pane, const std::string& chrom,
                            long long& s, long long& e, std::vector<float>& bins, float& mx) {
            float plotL = pane.x + 8, plotR = pane.getRight() - 8, plotW = std::max(1.f, plotR - plotL);
            long long span = e - s; if (span < 1) return;
            long long L = arcSize(chrom);
            double factor = (sy > 0) ? 0.8 : 1.25;
            long long newSpan = std::max(50LL, (long long)llround(span * factor));
            if (L > 0) newSpan = std::min(newSpan, L);
            float frac = ofClamp((x - plotL) / plotW, 0.f, 1.f);
            long long anchor = s + (long long)(frac * span);
            long long ns = anchor - (long long)(frac * newSpan), ne = ns + newSpan;
            if (ns < 0) { ns = 0; ne = newSpan; }
            if (L > 0 && ne > L) { ne = L; ns = std::max(0LL, L - newSpan); }
            s = ns; e = ne;
            bins = sumBinned(chrom, s, e, gwBinN);
            mx = 1.f; for (float v : bins) mx = std::max(mx, v);
        };
        if      (gwTopPane.inside(x, y)) { zoomPane(gwTopPane, gwTopChrom, gwTopS, gwTopE, gwTopBins, gwTopMax); requestGwPhase(); }
        else if (gwBotPane.inside(x, y)) { zoomPane(gwBotPane, gwBotChrom, gwBotS, gwBotE, gwBotBins, gwBotMax); requestGwPhase(); }
        return;
    }
    if (mode == Mode::Focus) {
        long long L = arcSize(focusOrigin.chrom), span = viewEnd - viewStart;
        double factor = (sy > 0) ? 0.8 : 1.25;
        long long newSpan = std::max(200LL, std::min((long long)llround(span*factor), L));
        float frac = ofClamp((angleAt(x, y) - fa0)/(fa1 - fa0), 0.f, 1.f);
        long long anchor = viewStart + (long long)(frac*span);
        long long ns = anchor - (long long)(frac*newSpan), ne = ns + newSpan;
        if (ns < 0) { ns = 0; ne = newSpan; }
        if (ne > L) { ne = L; ns = L - newSpan; }
        viewStart = ns; viewEnd = ne; rebinFocus();
    } else if (y >= layoutTop()) {
        int visible = (int)((ofGetHeight() - (layoutTop() + 34))/rowH);
        int maxS = std::max(0, (int)tableRows.size() - visible);
        tableScroll = (int)ofClamp(tableScroll - (int)(sy*3), 0, maxS);
    }
}

void ofApp::keyPressed(int key) {
    if (searchFocused) {                                 // typing into the search box
        if (key == OF_KEY_BACKSPACE) { if (!searchQuery.empty()) searchQuery.pop_back(); applySearch(); }
        else if (key == OF_KEY_RETURN || key == OF_KEY_ESC) searchFocused = false;
        else if (key >= 32 && key < 127) { searchQuery += (char)key; applySearch(); }
        return;
    }
    if (gwOpen) { if (key == OF_KEY_ESC) closeGW(); return; }             // GW overlay eats Esc first
    if (key == OF_KEY_ESC && mode == Mode::Focus) mode = Mode::Overview;  // pop zoom, keep trigger/arrows
    else if (key == OF_KEY_BACKSPACE && tableMode == TableMode::ByTrigger) backToChromTable();
}

void ofApp::windowResized(int w, int h) { recomputeLayout(); }
