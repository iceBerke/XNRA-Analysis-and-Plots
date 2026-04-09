import xml.etree.ElementTree as ET

file_path = r"C:\Users\berke\Desktop\Giuseppe\urgent_check\H392-15 550 °C NA 10 mg 166°.xnra"

IDF = "http://idf.schemas.itn.pt"
SIM = "http://www.simnra.com/simnra"
ns = {"idf": IDF, "simnra": SIM}

tree = ET.parse(file_path)
root = tree.getroot()

# Build parent pointers (read-only)
parent = {c: p for p in root.iter() for c in p}

def local(tag: str) -> str:
    # "{ns}name" -> "name"
    return tag.split("}", 1)[-1] if "}" in tag else tag

def xpath_of(elem) -> str:
    # Build an absolute-ish path with tag names and sibling indices
    parts = []
    cur = elem
    while cur is not None:
        p = parent.get(cur)
        if p is None:
            parts.append(local(cur.tag))
            break

        # index among same-tag siblings for disambiguation
        same = [c for c in list(p) if c.tag == cur.tag]
        if len(same) > 1:
            idx = same.index(cur) + 1  # 1-based
            parts.append(f"{local(cur.tag)}[{idx}]")
        else:
            parts.append(local(cur.tag))

        cur = p

    return "/" + "/".join(reversed(parts))

# ---- Your strategy, but with reporting ----
candidates = []

# Path 1
candidates.append((
    "path 1: .//simnra:calculateddata/idf:simpledata",
    root.find(".//simnra:calculateddata/idf:simpledata", ns)
))

# Path 2
candidates.append((
    "path 2: .//simnra:simulateddata/idf:simpledata",
    root.find(".//simnra:simulateddata/idf:simpledata", ns)
))

# Path 3
calc_spectrum = root.find('.//simnra:spectrum[@type="calculated"]', ns)
sd3 = None
if calc_spectrum is not None:
    sd3 = calc_spectrum.find(".//idf:simpledata", ns)
candidates.append((
    'path 3: .//simnra:spectrum[@type="calculated"]//idf:simpledata',
    sd3
))

# Path 4 (heuristic)
sd4 = None
hit_parent = None
for elem in root.iter():
    t = elem.tag.lower()
    if "calc" in t or "simul" in t:
        simdata = elem.find(".//idf:simpledata", ns)
        if simdata is not None:
            sd4 = simdata
            hit_parent = elem
            break
candidates.append((
    'path 4: first elem with "calc"/"simul" in tag, then .//idf:simpledata under it',
    sd4
))

# Pick the first match (same behavior as your code)
chosen_label = None
chosen = None
for label, elem in candidates:
    if elem is not None:
        chosen_label = label
        chosen = elem
        break

if chosen is None:
    print("No simulated_data found with the current candidate paths.")
else:
    print("MATCHED VIA:", chosen_label)
    print("simpledata path:", xpath_of(chosen))

    # Helpful: show where x/y live too
    x = chosen.find("idf:x", ns)
    y = chosen.find("idf:y", ns)
    if x is not None:
        print("x path:", xpath_of(x))
    if y is not None:
        print("y path:", xpath_of(y))

    # For path 4, also show which parent triggered the match
    if chosen_label.startswith("path 4") and hit_parent is not None:
        print("heuristic parent path:", xpath_of(hit_parent))
