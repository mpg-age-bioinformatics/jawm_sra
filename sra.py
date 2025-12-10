import jawm
import os
import gzip
import requests
from collections import Counter
import xml.etree.ElementTree as ET
import csv
from io import StringIO
from typing import Optional, Set, List
# from bs4 import BeautifulSoup

AGEPY_IMAGE="mpgagebioinformatics/agepy:9acb5de"

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _headers(email: Optional[str] = None) -> dict:
    return {
        "User-Agent": f"python-script (contact: {email or 'no-email-provided'})"
    }


def _get_runinfo_csv_for_bioproject(
    bioproject: str,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) -> str:
    """
    Use E-utilities to get the SRA RunInfo CSV for a BioProject (PRJNA...).
    Returns the raw CSV text.
    """
    # 1) ESEARCH to get WebEnv + QueryKey for the SRA records
    esearch_params = {
        "db": "sra",
        "term": bioproject,
        "usehistory": "y",
        "retmax": "0",
        "retmode": "xml",
    }
    if api_key:
        esearch_params["api_key"] = api_key

    r = requests.get(
        EUTILS_BASE + "esearch.fcgi",
        params=esearch_params,
        headers=_headers(email),
    )
    r.raise_for_status()

    root = ET.fromstring(r.text)
    webenv = root.findtext(".//WebEnv")
    query_key = root.findtext(".//QueryKey")

    if not webenv or not query_key:
        raise RuntimeError(
            f"Could not get WebEnv/QueryKey for {bioproject}. "
            f"Response (first 300 chars): {r.text[:300]!r}"
        )

    # 2) EFETCH in "runinfo" mode → CSV
    efetch_params = {
        "db": "sra",
        "query_key": query_key,
        "WebEnv": webenv,
        "rettype": "runinfo",
        "retmode": "text",
    }
    if api_key:
        efetch_params["api_key"] = api_key

    r2 = requests.get(
        EUTILS_BASE + "efetch.fcgi",
        params=efetch_params,
        headers=_headers(email),
    )
    r2.raise_for_status()
    csv_text = r2.text.strip()

    if not csv_text.startswith("Run,"):
        raise RuntimeError(
            f"RunInfo CSV for {bioproject} looks strange. "
            f"First 300 chars: {csv_text[:300]!r}"
        )

    return csv_text


def get_bioproject_organisms_from_sra(
    bioproject: str,
    email: Optional[str] = None,
    api_key: Optional[str] = None,
) :
    """
    Given a BioProject accession (e.g. 'PRJNA532357'),
    return a list of unique organism names seen in SRA RunInfo.

    Uses only standard library + 'requests', no Biopython, no pandas.
    """
    csv_text = _get_runinfo_csv_for_bioproject(
        bioproject, email=email, api_key=api_key
    )

    organisms: Set[str] = set()

    # Parse CSV using the standard library
    reader = csv.DictReader(StringIO(csv_text))
    # Common column names in RunInfo: 'ScientificName', sometimes 'Organism'
    for row in reader:
        name = (
            row.get("ScientificName")
            or row.get("Organism")
            or row.get("sample_name")
        )
        if name:
            name = name.strip()
            if name:
                organisms.add(name)

    if not organisms:
        raise ValueError(
            f"No organism names found in RunInfo for {bioproject} "
            f"(columns: {reader.fieldnames})"
        )

    # Return as a sorted list to have deterministic order
    return sorted(organisms)

# def get_geo_page_organism(accession: str) -> str:
#     """
#     Given a GEO accession like 'GSE129642', fetch the NCBI GEO page
#     and extract the 'Organism' entry (e.g. 'Caenorhabditis elegans').
#     """
#     accession = accession.strip()
#     url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"

#     resp = requests.get(url)
#     resp.raise_for_status()

#     soup = BeautifulSoup(resp.text, "html.parser")

#     # Find the row in the main info table where first cell == "Organism"
#     for row in soup.find_all("tr"):
#         cells = row.find_all(["td", "th"])
#         if len(cells) >= 2 and cells[0].get_text(strip=True) == "Organism":
#             organism = cells[1].get_text(strip=True)
#             return organism

#     raise ValueError(f"Could not find 'Organism' entry on page for {accession!r}")

def get_unique_sample_organism(accession: str):
    organisms = set()

    # if accession.startswith("GSE"):

    #     organism=get_geo_page_organism(accession)
    #     organisms.add( organism )

    # if accession.startswith("PRJ"):
    organisms = get_bioproject_organisms_from_sra(accession)

    if not organisms:
        raise ValueError("No organism found.")

    if len(organisms) > 1:
        raise ValueError(f"Multiple organisms found: {organisms}")

    return organisms.pop().lower().replace(" ", "_")


def get_nconcatenations( tsv ) :
    first_col = []
    with open(tsv, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if parts:  # skip empty lines
                first_col.append(parts[0])
                if parts[4] == "PAIRED" :
                    first_col.append(parts[0])

    # Count occurrences
    counts = Counter(first_col)

    num_duplicated_items = sum(1 for c in counts.values() if c > 1)

    return num_duplicated_items

read_acc=jawm.Process( 
    name="read_acc",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"], f"{p.var['acc']}.samples.xlsx"  ) ) ,
    script="""#!/usr/bin/env python3
import AGEpy as age
import os
import requests
import pandas as pd
import re
import unicodedata

def url_exists(url: str, timeout: int = 5) -> bool:
    try:
        # Try HEAD first (lighter)
        resp = requests.head(url, allow_redirects=True, timeout=timeout)
        if resp.status_code >= 400:
            # Some servers don't support HEAD properly; fall back to GET
            resp = requests.get(url, stream=True, timeout=timeout)
        return resp.ok
    except requests.RequestException:
        return False

def load_table(file_path):
    _, ext = os.path.splitext(file_path)
    ext = ext.lower()

    if ext in [".xls", ".xlsx"]:
        # Excel file
        df=pd.read_excel(file_path)
    
    elif ext in [".tsv", ".csv"]:

        # TSV file (tab-separated)
        print("Reading file as tab separated file")

        try:
            
            df=pd.read_csv( file_path, sep="\\t" )

        except:
            
            raise ValueError(f"Could not read tab separated file")

    else:
        
        raise ValueError(f"Unsupported file type: {ext}")

    if len( df.columns.tolist() ) < 2 :
        raise ValueError(f"Imported table had less than 2 columns.")

    if len( df.columns.tolist() ) > 2 :

        cols = list(df.columns)
        cols_=[ s for s in cols if s == "Experiment" ]+[ s for s in cols if s in ["group", "groups"] ]

        if len(cols_) == 2 :
            df=df[cols_]

        else:

            exp_idx = cols.index("Experiment")

            # Check that a column exists after Experiment
            if exp_idx + 1 >= len(cols):
                raise ValueError("Column 'Experiment' exists but there is no column after it.")

            # Select Experiment and the following column
            df=df[ cols[exp_idx : exp_idx + 2] ]
    
    df.columns=["sample","group"]

    return df



def process_groups(groups):

    if os.path.isfile( groups ) :
        groups=load_table( file_path )

    else:

        # reading from a variable of the form "sample1:group;sample1:group;.."
        groups=groups.replace("\\x01", "" )
        groups=groups.strip().strip(";")
        groups = re.sub(r'\s*([;])\s*', r'\\1', groups )
        groups=[ s.split(";") for s in groups.split("\\n") ]
        groups=pd.DataFrame(groups, columns=["sample","group"] )

    for c in groups.columns.tolist():
        groups[c]=groups[c].apply(lambda x: x.strip() )
    groups=groups.drop_duplicates(subset=["sample"])
    groups["group"] = groups["group"].apply(lambda x: age.safe_filename(x) )
    return groups

if url_exists("{{url_datasheet}}"):
    
    # reading pre-ready datasheets from an url
    samples_df = pd.read_csv("{{url_datasheet}}", sep="\\t")

elif ( "{{groups}}" == "" ) and "{{acc}}".startswith("GSE") :

    # generating samples sheet from a geo accession number
    samples_df, groups_df, runinfo_df= age.geo_to_sra_samples( "{{acc}}", email="{{email}}", api_key="{{api_key}}" )

elif ( "{{groups}}" == "" ) and "{{acc}}".startswith("PRJ") :

    raise ValueError("BioProject accessions must be supplemented with group information for each sample.")

# elif "{{acc}}".startswith("GSE") :

#     groups=process_groups( "{{groups}}" )
#     gsms=groups["sample"].tolist()
#     samples_df=age.gsms_to_sra_samples( gsms,  email="{{email}}", api_key="{{api_key}}" )
#     samples_df=samples_df.drop( [ "group" ] , axis=1 )
#     samples_df=pd.merge( groups, samples_df, on=["sample"], how="inner" )
#     cols = list(samples_df.columns)
#     col = cols.pop(7)
#     cols.insert(5, col)
#     samples_df = samples_df[cols]

elif "{{acc}}".startswith("PRJ") :

    groups=process_groups("{{groups}}")
    samples_df=age.fetch_sra_metadata_table_for_bioproject( "{{acc}}" )[[ "Run", "Experiment", "LibraryLayout" ]]
    samples_df=pd.merge( groups, samples_df, left_on=[ "sample" ], right_on=[ "Experiment" ], how="inner" )

if not os.path.isdir("{{raw_data}}" ) :
    os.makedirs( "{{raw_data}}" )

outfile_xlsx=os.path.join( "{{raw_data}}" , "{{acc}}"+".samples.xlsx" )
outfile_tsv=os.path.join( "{{raw_data}}" , "{{acc}}"+".samples.tsv" )

samples_df.to_excel( outfile_xlsx, index=False )
samples_df.to_csv( outfile_tsv, sep="\\t", index=False )
""",
    var={
        "email":"",
        "api_key":""
        },
    desc={
        "email" : "email passed to NCBI E-utilities",
        "api_key" : "NCBI API key for higher request limits",
        "raw_data" : "Downloads folder.",
        "acc" : "geo accession or sra bioproject",
        "url_datasheet" : "",
        "groups" : "Instead of giving in an accession number or a preset sample sheet you can also\
            supply a table in string form '<geo_sample>;<group_name>\n<geo_sample>;<group_name>\n..' \
            otherwise leave it as ''"
    },
    container=AGEPY_IMAGE
)

prefetch=jawm.Process( 
    name="prefetch",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"], "sra", f'{p.var["sraid"]}.touch'  ) ) ,
    script="""#!/bin/bash
cd {{raw_data}}
mkdir sra
cd sra
prefetch --max-size u -p {{sraid}}
while [ $? -ne 0 ]; do 
   prefetch --max-size u -p {{sraid}}
done
""",
    var={},
    desc={
        "raw_data": "Downloads folder.",
        "sraid":"sra_id", 
    },
    manager_slurm={"--mem":"20GB", "-t":"1:00:00", "-c":"8" },
    container="mpgagebioinformatics/sra:3.2.1"
)

fastq_dump=jawm.Process( 
    name="fastq_dump",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"], "sra", f'{p.var["sraid"]}.touch'  ) ) ,
    script="""#!/bin/bash
set -e
cd {{raw_data}}/sra

shopt -s nullglob

for f in {{sraid}}*lite ; do
    fastq-dump --split-3 "./${f}"
    rm -rf "${f}"
done

if [ -d {{sraid}} ] ; then
    fastq-dump --split-3 "./{{sraid}}"
    rm -rf "{{sraid}}"
fi

if [ -f {{sraid}}.fastq ] ; then pigz {{sraid}}.fastq ; fi
if [ -f {{sraid}}_1.fastq ] ; then pigz {{sraid}}_1.fastq ; fi
if [ -f {{sraid}}_2.fastq ] ; then pigz {{sraid}}_2.fastq ; fi

touch {{sraid}}.touch
rm -rf {{sraid}}
""",
    var={},
    desc={
        "raw_data": "Downloads folder.",
        "sraid":"sra_id", 
    },
    manager_slurm={"--mem":"20GB", "-t":"1:00:00", "-c":"8" },
    container="mpgagebioinformatics/sra:3.2.1"
)


relabel_sra=jawm.Process( 
    name="relabel_sra",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"], "sra", "relabel_sra.touch"  ) ) ,
    script="""#!/usr/bin/env python3
import pandas as pd
import os
import gzip
import shutil
from pathlib import Path
import re
import unicodedata
import AGEpy as age
from concurrent.futures import ProcessPoolExecutor, as_completed


def concat_gz_files( files, output, concat="{{concat}}" ):

    if concat == "text": 
        with gzip.open(output, "wt") as outfile :
            for fname in files:
                with gzip.open(fname, "rt") as infile :
                    for line in infile:
                        outfile.write(line)

    elif concat == "stream" :
        with open(output, "wb") as outfile :
            for fname in files:
                with open(fname, "rb") as infile :
                    shutil.copyfileobj(infile, outfile, length=1024 * 1024)



def add_rep_suffix(groups):
    counts = {}
    output = []
    for g in groups:
        counts[g] = counts.get(g, 0) + 1
        output.append(f"rep_{counts[g]}")
    return output


raw_data = os.path.abspath("{{raw_data}}")

input_xlsx = os.path.join(raw_data, "{{acc}}" + ".samples.xlsx")

os.chdir(f"{{raw_data}}/sra")

df = pd.read_excel(input_xlsx)
df = (
    df.groupby("sample", as_index=False)
      .agg(lambda s: ",".join(sorted(set(s.dropna().astype(str)))))
)

if ( "{{groups}}" != "" ) and ( not os.path.isfile( "{{groups}}" )  ):
    # reading from a variable of the form "sample1:group;sample1:group;.."
    groups = "{{groups}}".replace("\\x01", "")
    groups = groups.strip().strip(";")
    groups = re.sub(r'\s*([;])\s*', r'\\1', groups)
    groups = [s.split(";") for s in groups.split("\\n")]
    groups = pd.DataFrame(groups, columns=["sample", "group"])
    groups = groups.drop(["group"], axis=1)
    groups["sample"] = groups["sample"].apply(lambda x: x.strip())
    df["sample"] = df["sample"].apply(lambda x: x.strip())
    df = pd.merge(groups, df, on=["sample"], how="inner")

df["group"] = df["group"].apply(lambda x: age.safe_filename(x))
df["rep"] = add_rep_suffix(df["group"])

# Collect concatenation jobs to run in parallel
futures = []

with ProcessPoolExecutor(max_workers=int({{parallel}})) as executor:
    for experiment in df["Experiment"].tolist():
        runs = df.loc[df["Experiment"] == experiment, "Run"].iloc[0].split(",")
        layout = df.loc[df["Experiment"] == experiment, "LibraryLayout"].iloc[0]
        rep = df.loc[df["Experiment"] == experiment, "rep"].iloc[0]
        group = df.loc[df["Experiment"] == experiment, "group"].iloc[0]

        # Single run: just rename synchronously
        if len(runs) == 1:
            run = runs[0]

            if layout == "PAIRED":
                old_name_1 = f"{run}_1.fastq.gz"
                new_name_1 = os.path.join(raw_data, f"{group}.{rep}.read_1.fastq.gz")
                print(f"Renaming {old_name_1} as {new_name_1}")
                os.rename(old_name_1, new_name_1)

                old_name_2 = f"{run}_2.fastq.gz"
                new_name_2 = os.path.join(raw_data, f"{group}.{rep}.read_2.fastq.gz")
                print(f"Renaming {old_name_2} as {new_name_2}")
                os.rename(old_name_2, new_name_2)

            else:
                old_name = f"{run}.fastq.gz"
                new_name = os.path.join(raw_data, f"{group}.{rep}.read_1.fastq.gz")
                print(f"Renaming {old_name} as {new_name}")
                os.rename(old_name, new_name)

        # Multiple runs: schedule concatenations in the process pool
        else:
            if layout == "PAIRED":
                new_name_1 = os.path.join(raw_data, f"{group}.{rep}.read_1.fastq.gz")
                new_name_2 = os.path.join(raw_data, f"{group}.{rep}.read_2.fastq.gz")

                files_1 = [f"{s}_1.fastq.gz" for s in runs]
                files_2 = [f"{s}_2.fastq.gz" for s in runs]

                files_1_str = ", ".join(files_1)
                files_2_str = ", ".join(files_2)

                print(f"Concatenating {files_1_str} as {new_name_1}")
                futures.append(executor.submit(concat_gz_files, files_1, new_name_1))

                print(f"Concatenating {files_2_str} as {new_name_2}")
                futures.append(executor.submit(concat_gz_files, files_2, new_name_2))

            else:
                new_name = os.path.join(raw_data, f"{group}.{rep}.read_1.fastq.gz")
                files = [f"{s}.fastq.gz" for s in runs]
                files_str = ", ".join(files)

                print(f"Concatenating {files_str} as {new_name}")
                futures.append(executor.submit(concat_gz_files, files, new_name))

    # Wait for all concatenations to complete and propagate any exceptions
    for fut in as_completed(futures):
        fut.result()

# Only touch the file once everything is finished
Path(os.path.join(raw_data, "sra", "relabel_sra.touch")).touch()
""",
    var={
        "parallel": 4,
        "concat":"stream"
    },
    desc={
        "concat": "'stream or text for how you want to concat fastq.gz files when multiple files per sample exist. Default='stream' ",
        "raw_data": "Downloads folder.",
        "acc":"geo accession",
        "groups":"Instead of giving in an accession number or a preset sample sheet you can also\
            supply a table in string form '<geo_sample>;<group_name>\n<geo_sample>;<group_name>\n..' \
            otherwise leave it as ''"    },
    manager_slurm={"--mem":"20GB", "-t":"1:00:00", "-c":"4" },
    container=AGEPY_IMAGE
)

# sorted.fastq.gz
test_unpigz=jawm.Process( 
    name="test_unpigz",
    when=lambda p: not os.path.isfile( os.path.join( p.var["raw_data"], f'{p.var["sraid"]}.fastq'  ) ) ,
    script="""#!/bin/bash
cd {{raw_data}}
unpigz {{sraid}}.fastq.gz
""",
    var={},
    desc={
        "raw_data": "Downloads folder.",
        "sraid":"sra_id", 
    },
    container="mpgagebioinformatics/sra:3.2.1"
)


if __name__ == "__main__":
    import sys
    from jawm.utils import workflow, load_modules, get_image

    # pre-download images to avoid parallel 
    # processes concurring to download the same image
    images=get_image( [ read_geo.container, prefetch.container ] )

    # parse arguments
    workflows, var, args, unknown_args=jawm.utils.parse_arguments( ["main","sra","test", "geo"] )

    if workflow( ["geo"], workflows ) :
        import pandas as pd

        # read the respective geo accession and create samples sheet
        read_geo.execute()
        jawm.Process.wait( read_geo.hash )

        input_xlsx=os.path.join( read_geo.var["raw_data"], f"{read_geo.var['acc']}.samples.xlsx"  )

        df=pd.read_excel( input_xlsx )
        sra_ids=df["Run"].tolist()

    else:

        sra_ids=[]

    # usage: 
    if workflow( ["main","sra","test"], workflows ) :  

        if not sra_ids : 
            sra_ids =  prefetch.var["sraid"].split(",")  

        fastq_dump_jobs=[]
        for sraid in sra_ids : 
            
            # download serially to prevent too much io and server query
            prefetch_=prefetch.clone()
            prefetch_.var["sraid"]=sraid
            prefetch_.execute()
            jawm.Process.wait( prefetch_.hash )

            # send extraction jobs to background
            fastq_dump_=fastq_dump.clone()
            fastq_dump_.var["sraid"]=sraid
            fastq_dump_.execute(  )

            fastq_dump_jobs.append( fastq_dump_.hash )
        
        jawm.Process.wait( fastq_dump_jobs )

    if workflow( ["geo"], workflows ) and workflow( ["sra"], workflows ) :

        relabel_sra.execute( )

        jawm.Process.wait( relabel_sra.hash )

    if workflow( "test", workflows ) :

        # for the test workflow we also do something more (just for sra)
        test_unpigz.execute()
        jawm.Process.wait( )
        print("Test completed.")

    sys.exit(0)
