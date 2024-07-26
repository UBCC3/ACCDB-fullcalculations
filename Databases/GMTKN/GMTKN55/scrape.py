from lxml import html
import requests
import csv

#URL_BASE = "http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/GMTKN55"
URL_BASE = "file://./GMTKN55"

indent = 0
def ipush():
    global indent
    indent += 1
def ipop():
    global indent
    indent -= 1
def iprint(*args):
    global indent
    print(("  "*indent)+args[0], *args[1:])

def write_csv(file, **kwargs):
    """
        Write a CSV file with headers matching the kwargs, and items matching the iterable arguments.
    """
    with open(file, 'w', newline='') as output_file:
        csv_data = csv.writer(output_file, dialect='excel')
        
        csv_data.writerow([*kwargs.keys()])
        for key in kwargs.keys(): kwargs[key] = kwargs[key].__iter__()
        
        n_keys = len(kwargs.keys())
        n_stop = 0
        while True:
            row = []
            for arg in kwargs.values():
                try: row.append(arg.__next__())
                except StopIteration: n_stop += 1
            if n_stop == n_keys: return
            if n_stop > 0: raise RuntimeError('Arguments to writeCsv have different lengths!')
            csv_data.writerow(row)

def fetch(url):
    if url[0:7].lower() == "file://":
        # Ok, not technically a URL anymore...
        with open(url[7:], 'r') as f:
            return f.read()
    else:
        page = requests.get(url)
        if page.status_code != 200: raise BaseException(page.status_code)
        return page.content

def get_rows(url):
    ipush()
    iprint("Scrape", url)

    tree = html.fromstring(fetch(url))
    tables = tree.findall('.//table')
    if len(tables) != 1: raise BaseException("Wrong number of tables " + str(len(tables)))
    rows = tables[0].xpath('tr')
    headers = [s.strip() for s in rows[0].xpath('th/text()')]
    rows = [[s.strip() for s in row.xpath('td/text()')] for row in rows[1:]]
    ipop()
    return (headers, rows)

def do_scrape(ds, dft):
    ipush()
    iprint("Fetching DFT", dft, "for dataset", ds)

    (headers, rows) = get_rows('{}/results/{}/{}/result.html'.format(URL_BASE, ds, dft))
    dft_names = headers[headers.index("Ref.") + 1:]

    dft_names[0] = dft
    dft_names[1:] = [dft_names[0]+"-"+n for n in dft_names[1:]]

    dp_map = {}
    for row in rows:
        if not len(row): continue
        if not row[0].isdigit(): break

        n = int(row[0])
        dp_name = "{}_{}".format(ds, n)

        # All values are expressed as deviation from reference
        i = -(len(dft_names) + 1)
        ref = float(row[i])
        for dftn in dft_names:
            i += 1
            if not row[i]:
                iprint("WARNING: Missing datapoint for point", dp_name, "with DFT", dftn)
                continue
            try:
                dp_map[(dp_name, dftn)] = ref+float(row[i])
            except ValueError:
                iprint("WARNING: Malformed datapoin", dp_name, "with DFT", dftn, "=", row[i])
    ipop()
    return dp_map

def do_scrape_reference(ds):
    ipush()
    iprint("Fetching reference for dataset", ds)
    dp_map = {}
    (headers, rows) = get_rows('{}/{}ref.html'.format(URL_BASE, ds))
    for row in rows:
        if not len(row): continue
        if not row[0].isdigit(): break

        n = int(row[0])
        dp_name = "{}_{}".format(ds, n)

        # All values are expressed as deviation from reference
        dp_map[dp_name] = float(row[-1])
    ipop()
    return dp_map

datasets = ["W4-11", "G21EA", "G21IP", "DIPCS10", "PA26", "SIE4x4", "ALKBDE10", "YBDE18", "AL2X6", "HEAVYSB11", "NBPRC", "ALK8", "RC21", "G2RC", "BH76RC", "FH51", "TAUT15", "DC13", "MB16-43", "DARC", "RSE43", "BSR36", "CDIE20", "ISO34", "ISOL24", "C60ISO", "PArel", "BH76", "BHPERI", "BHDIV10", "INV24", "BHROT27", "PX13", "WCPT18", "RG18", "ADIM6", "S22", "S66", "HEAVY28", "WATER27", "CARBHB12", "PNICO23", "HAL59", "AHB21", "CHB6", "IL16", "IDISP", "ICONF", "ACONF", "Amino20x4", "PCONF21", "MCONF", "SCONF", "UPU23", "BUT14DIOL"]
dfts = ["PBE", "PBEhPBE", "revPBE", "RPBE", "PW91", "BLYP", "BP86", "BPBE", "OPBE", "OLYP", "XLYP", "mPWLYP", "PW91P86", "mPWPW91", "rPW86PBE", "B97-D3", "HCTH", "N12", "VV10", "PKZB", "TPSS", "revTPSS", "SCAN", "tHCTH", "M06L", "M11L", "MN12L", "MN15L", "B3LYP", "B3LYP-NL", "B3PW91", "B3P86", "BHLYP", "B1P86", "B1LYP", "B1B95", "MPW1B95", "PW6B95", "MPWB1K", "mPW1LYP", "mPW1PW91", "PW1PW", "MPW1KCIS", "MPWKCIS1K", "PBE0", "PBEh1PBE", "PBE1KCIS", "X3LYP", "O3LYP", "B97-1", "B97-2", "B98", "HISS", "HSE03", "HSE06", "TPSSh", "revTPSSh", "TPSS0", "revTPSS0", "TPSS1KCIS", "BMK", "tHCTHhyb", "M05", "M052X", "M06", "M062X", "M08HX", "M11", "SOGGA11X", "N12SX", "MN12SX", "MN15", "LC-whPBE", "wB97X-D3", "wB97X-V", "APFD", "B2PLYP", "B2GPPLYP", "MPW2PLYP", "PWPB95", "DSD-BLYP", "DSD-PBEP86", "DSD-PBEB95", "r2SCAN-3c"]

dfts = sorted(dfts)

# Use with python3 scrape.py | parallel -j 10 --gnu curl 'http://www.thch.uni-bonn.de/tc.old/downloads/GMTKN/{}' --create-dirs -o {}
#for ds in datasets:
#    print('GMTKN55/{}ref.html'.format(ds))
#    for dft in dfts:
#        print('GMTKN55/results/{}/{}/result.html'.format(ds, dft))
#exit(0)

big_map = {}
for ds in datasets:
    iprint("Dataset", ds)
    smaller_map = {(dp, "Reference"): v for dp, v in do_scrape_reference(ds).items()}
    for dft in dfts:
        ipush()
        iprint("DFT", dft)
        smaller_map = {**smaller_map, **do_scrape(ds, dft)}
        ipop()
    #iprint("Smaller map is:", smaller_map)
    big_map = {**big_map, **smaller_map}

#iprint("Big map is:", big_map)
iprint("Writing dataset")

dps = map(lambda n: n[0], big_map.keys())
dps = sorted(set(dps), key=lambda n: ("_".join(n.split('_')[:-1]), int(n.split('_')[-1])))

ddfts = map(lambda n: n[1], big_map.keys())
ddfts = sorted(set(ddfts))

malformed_dps = []
def lookup(dp, dft):
    global malformed_dps
    # Some datapoints are malformed
    if (dp, dft) in big_map:
        return "{:.5f}".format(big_map[(dp, dft)])
    else:
        malformed_dps.append((dp, dft))
        return "??"

write_csv(
    "IndValues.csv",
    RefNames        = dps,
    DatasetRefNames = ["_".join(n.split('_')[:-1]) for n in dps],
    RefValues       = [lookup(n, "Reference") for n in dps],
    **{dft: [lookup(dp, dft) for dp in dps] for dft in ddfts}
)

print()
print("Malformed Datapoints:")
for (dp, dft) in malformed_dps: print(dp, dft)