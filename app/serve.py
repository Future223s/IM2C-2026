"""Start the Peak Foliage Explorer on http://localhost:8000 (standard library only, any Python 3).

    python app/serve.py            # then open http://localhost:8000
    python app/serve.py --port 8080
"""
import argparse
import functools
import http.server
import webbrowser
from pathlib import Path

HERE = Path(__file__).resolve().parent

ap = argparse.ArgumentParser()
ap.add_argument("--port", type=int, default=8000)
ap.add_argument("--no-browser", action="store_true")
args = ap.parse_args()

Handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(HERE))
with http.server.ThreadingHTTPServer(("127.0.0.1", args.port), Handler) as httpd:
    url = f"http://localhost:{args.port}"
    print(f"Peak Foliage Explorer running at {url}  (Ctrl+C to stop)")
    if not args.no_browser:
        webbrowser.open(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
