"""Build stamp: which code, seed, constants and inputs made the output in $YC_WORK.

  python scripts/build_stamp.py start [FIRST_STAGE]   write (stage 1) or extend (restart from a later stage) $WORK/build.json
  python scripts/build_stamp.py check                 list result files in $WORK that are older than the build, and
                                                      write $WORK/validation/MANIFEST.json (file -> build id, age)
  python scripts/build_stamp.py id                    print the build id

build.json: build_id, commit, dirty (uncommitted changes in pipeline/, scripts/, validation/), seed, start time, host,
workers, every module-level constant of pipeline/*.py (literal values; other expressions as source text) and size and
modification time of the inputs. Why: with several builds a day, numbers from an earlier build reached the documents
twice (2026-11-12); a result file without the current build id, or older than the build, is not a result of it.
"""
import ast, glob, json, os, socket, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(BUNDLE, "pipeline"))
import config as C  # noqa: E402

F = f"{C.WORK}/build.json"


def _git(*a):
    try:
        return subprocess.run(["git", "-C", BUNDLE, *a], capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return ""


def constants():
    """module-level UPPER_CASE assignments of every pipeline module, without importing it"""
    out = {}
    for p in sorted(glob.glob(f"{BUNDLE}/pipeline/*.py")):
        src = open(p).read()
        try:
            tree = ast.parse(src)
        except SyntaxError:
            continue
        d = {}
        for node in tree.body:
            if isinstance(node, ast.Assign):
                names = [t.id for t in node.targets if isinstance(t, ast.Name)]
                if isinstance(node.targets[0], ast.Tuple):
                    names = [e.id for e in node.targets[0].elts if isinstance(e, ast.Name)]
            elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None:
                names = [node.target.id]
            else:
                continue
            names = [n for n in names if n.isupper() or (n.replace("_", "").isupper() and not n.startswith("_"))]
            if not names:
                continue
            try:
                v = ast.literal_eval(node.value)
                json.dumps(v)
            except Exception:
                v = "= " + ast.get_source_segment(src, node.value)[:300]
            d[", ".join(names)] = v
        if d:
            out[os.path.basename(p)[:-3]] = d
    return out


def inputs():
    """size and modification time of the inputs (data files, library size, large downloads)"""
    out = {}
    for p in sorted(glob.glob(f"{C.DATA}/*")) + [C.AFDB_TAR, C.POOL_SCORES, f"{C.WORK}/rcsb_cache.json",
                                                 f"{C.ROOT}/cell_meshes.npz", f"{C.ROOT}/organelle_stats.json"]:
        if os.path.isfile(p):
            st = os.stat(p)
            out[os.path.relpath(p, BUNDLE) if p.startswith(BUNDLE) else p] = [st.st_size, time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime))]
    lib = glob.glob(f"{C.LIB}/*.npz")
    out["pipeline/lib (npz files, newest)"] = [len(lib), time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(max(map(os.path.getmtime, lib)))) if lib else None]
    return out


def start(first=1):
    now = time.strftime("%Y-%m-%dT%H:%M:%S")
    sha = _git("rev-parse", "HEAD")
    dirty = bool(_git("status", "--porcelain", "--untracked-files=no", "--", "pipeline", "scripts", "validation"))
    if first > 1 and os.path.exists(F):
        b = json.load(open(F))
        b.setdefault("restarts", []).append(dict(time=now, first_stage=first, commit=sha, dirty=dirty, seed=C.SEED))
        if sha != b.get("commit") or dirty or C.SEED != b.get("seed"):
            b["mixed"] = True                 # later stages ran with other code or another seed than the first ones
    else:
        b = dict(build_id=f"{time.strftime('%Y%m%d-%H%M%S')}-{sha[:7] or 'nogit'}{'-dirty' if dirty else ''}-s{C.SEED}",
                 commit=sha, dirty=dirty, seed=C.SEED, started=now, t0=time.time(), host=socket.gethostname(), nproc=C.NPROC,
                 first_stage=first, constants=constants(), inputs=inputs())
    json.dump(b, open(F, "w"), indent=1)
    print("build", b["build_id"], "(mixed)" if b.get("mixed") else "")


def build_id():
    return json.load(open(F)).get("build_id") if os.path.exists(F) else None


def check():
    if not os.path.exists(F):
        print("no build.json in", C.WORK); return 1
    b = json.load(open(F)); t0 = b["t0"]
    files = sorted(glob.glob(f"{C.WORK}/validation/*.json") + glob.glob(f"{C.WORK}/validation/latest/*.json")
                   + glob.glob(f"{C.RESULTS}/*.json") + glob.glob(f"{C.WORK}/condensates/*.json"))
    man, stale = {}, []
    for p in files:
        if os.path.basename(p) == "MANIFEST.json":
            continue
        mt = os.path.getmtime(p); rel = os.path.relpath(p, C.WORK)
        ok = mt >= t0
        man[rel] = dict(modified=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(mt)), of_this_build=ok)
        if not ok:
            stale.append(rel)
    os.makedirs(f"{C.WORK}/validation", exist_ok=True)
    json.dump(dict(build_id=b["build_id"], mixed=b.get("mixed", False), checked=time.strftime("%Y-%m-%dT%H:%M:%S"), files=man),
              open(f"{C.WORK}/validation/MANIFEST.json", "w"), indent=1)
    print(f"build {b['build_id']}: {len(man) - len(stale)} result files of this build, {len(stale)} older than it")
    for s in stale:
        print("  STALE", man[s]["modified"], s)
    return 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "id"
    if cmd == "start":
        start(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
    elif cmd == "check":
        sys.exit(check())
    else:
        print(build_id())
