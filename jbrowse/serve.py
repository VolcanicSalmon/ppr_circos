#!/usr/bin/env python3
"""
Local JBrowse server for the trackplot 'JB' button.

- Serves this jbrowse/ dir statically (config.json, tracks/, and a JBrowse2 web build).
- Reverse-proxies the remote reference genomes (which send NO CORS) so JBrowse can stream
  them same-origin:  /ref/rh/*  -> SpudDB (RH_v3),  /ref/nb/* -> lifenglab (NBE_HZ).
  Range: headers are forwarded, so genomes are streamed (never downloaded whole), and the
  proxied responses get Access-Control-Allow-Origin so the browser allows them.

usage:  python3 serve.py [port]        (default 9000)
then open the URL the app's JB button produces, e.g. http://localhost:9000/?assembly=RH_v3&loc=...
"""
import sys, os, json, urllib.request, urllib.error
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer

PORT = int(sys.argv[1]) if len(sys.argv) > 1 else 9000
HERE = os.path.dirname(os.path.abspath(__file__))

# Built-in genome proxies. add_jbrowse.py appends user entries to proxies.json,
# which is merged over these at startup (so the helper never edits this file).
PROXIES = {
    "/ref/rh/":   "https://spuddb.uga.edu/jb2/",
    "/ref/nb/":   "http://lifenglab.hzau.edu.cn/Nicomics/Download/NbeHZ1/",
    "/ref/spim/": "https://solgenomics.net/ftp/genomes/Solanum_pimpinellifolium/LA2093/Spimp_LA2093_genome_v1.5/",
    "/ref/aras/": "https://jbrowse2.arabidopsis.org/",
}
_pj = os.path.join(HERE, "proxies.json")
if os.path.exists(_pj):
    try:
        with open(_pj) as _fh:
            PROXIES.update(json.load(_fh))
    except Exception as _e:
        print(f"warning: could not read proxies.json: {_e}")


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **k):
        super().__init__(*a, directory=HERE, **k)

    def _proxy_target(self):
        for prefix, base in PROXIES.items():
            if self.path.startswith(prefix):
                return base + self.path[len(prefix):]
        return None

    def _do_proxy(self, method):
        target = self._proxy_target()
        if target is None:
            return False
        headers = {}
        if "Range" in self.headers:
            headers["Range"] = self.headers["Range"]
        headers["User-Agent"] = "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15) trackplot-jb-proxy"
        req = urllib.request.Request(target, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                self.send_response(r.status)
                for h in ("Content-Type", "Content-Length", "Content-Range", "Accept-Ranges", "Last-Modified", "ETag"):
                    v = r.headers.get(h)
                    if v:
                        self.send_header(h, v)
                self.send_header("Access-Control-Allow-Origin", "*")
                self.send_header("Access-Control-Allow-Headers", "Range")
                self.end_headers()
                if method == "GET":
                    while True:
                        chunk = r.read(65536)
                        if not chunk:
                            break
                        self.wfile.write(chunk)
        except urllib.error.HTTPError as e:
            self.send_response(e.code); self.send_header("Access-Control-Allow-Origin", "*"); self.end_headers()
        except Exception as e:
            self.send_error(502, f"proxy error: {e}")
        return True

    def do_GET(self):
        if not self._do_proxy("GET"):
            super().do_GET()

    def do_HEAD(self):
        if not self._do_proxy("HEAD"):
            super().do_HEAD()

    def end_headers(self):
        # CORS for the static (tracks/) responses too, harmless
        if not self._proxy_target():
            self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()


if __name__ == "__main__":
    os.chdir(HERE)
    print(f"serving {HERE} at http://localhost:{PORT}  (proxying /ref/rh -> SpudDB, /ref/nb -> lifenglab)")
    ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
