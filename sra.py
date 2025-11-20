import jawm
import os


read_geo=jawm.Process( 
    name="read_geo",
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], f"{p.var['geoacc']}.samples.xlsx"  ) ) ,
    script="""#!/usr/bin/env python3
import AGEpy as age
import os
import requests
import pandas as pd

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

if url_exists("{{url_datasheet}}"):
    
    # reading pre-ready datasheets from an url
    samples_df = pd.read_csv("{{url_datasheet}}", sep="\\t")

elif "{{groups}}" == "" :

    # generating samples sheet from a geo accession number
    samples_df, groups_df, runinfo_df= age.geo_to_sra_samples( "{{geoacc}}" )

else :

    # reading from a variable of the form "sample1:group;sample1:group;.."
    groups=[ s.split(":") for s in "{{groups}}".replace(" ", "").split(";") ]
    groups=pd.DataFrame(groups, columns=["sample","group"] )
    gsms=groups["sample"].tolist()
    samples_df=age.gsms_to_sra_samples(gsms)
    samples_df=samples_df.drop( [ "group" ] , axis=1 )
    samples_df=pd.merge( groups, samples_df, on=["sample"], how="inner" )


if not os.path.isdir("{{downloads_folder}}" ) :
    os.makedirs( "{{downloads_folder}}" )

outfile_xlsx=os.path.join( "{{downloads_folder}}" , "{{geoacc}}"+".samples.xlsx" )
outfile_tsv=os.path.join( "{{downloads_folder}}" , "{{geoacc}}"+".samples.tsv" )

samples_df.to_excel( outfile_xlsx, index=False )
samples_df.to_csv( outfile_tsv, sep="\\t", index=False )
""",
    desc={
        "downloads_folder": "Downloads folder.",
        "geoacc":"geo accession",
        "url_datasheet":"",
        "groups":"Instead of giving in an accession number or a preset sample sheet you can also\
            supply a dictionary like '<geo_sample>:<group_name>;<geo_sample>:<group_name>;..' \
            otherwise leave it as ''"
    },
    container="mpgagebioinformatics/agepy:f92ee75"
)

prefetch=jawm.Process( 
    name="prefetch",
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], f'{p.var["sraid"]}.touch'  ) ) ,
    script="""#!/bin/bash
cd {{downloads_folder}}
mkdir sra
cd sra
prefetch --max-size u -p {{sraid}}
while [ $? -ne 0 ]; do 
   prefetch --max-size u -p {{sraid}}
done
""",
    desc={
        "downloads_folder": "Downloads folder.",
        "sraid":"sra_id", 
    },
    container="mpgagebioinformatics/sra:3.2.1"
)

fastq_dump=jawm.Process( 
    name="fastq_dump",
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], f'{p.var["sraid"]}.touch'  ) ) ,
    script="""#!/bin/bash
set -e
cd {{downloads_folder}}/sra

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
    desc={
        "downloads_folder": "Downloads folder.",
        "sraid":"sra_id", 
    },
    container="mpgagebioinformatics/sra:3.2.1"
)


relabel_geo=jawm.Process( 
    name="relabel_geo",
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], "relabel_geo.touch"  ) ) ,
    script="""#!/usr/bin/env python3
import pandas as pd
import os
import gzip
import shutil
from pathlib import Path

def concat_gz_files(files, output):
    with gzip.open(output, "wt") as outfile:
        for fname in files:
            with gzip.open(fname, "rt") as infile:
                for line in infile:
                    outfile.write(line)

def add_rep_suffix(groups):
    counts = {}
    output = []
    for g in groups:
        counts[g] = counts.get(g, 0) + 1
        output.append(f"rep_{counts[g]}")
    return output

os.chdir(f"{{downloads_folder}}/sra")

input_xlsx=os.path.join( "{{downloads_folder}}" , "{{geoacc}}"+".samples.xlsx" )

df=pd.read_excel( input_xlsx )
df = (
    df.groupby("sample", as_index=False)
      .agg(lambda s: ",".join(sorted(set(s.dropna().astype(str)))))
)

if "{{groups}}" != ""  :

    # reading from a variable of the form "sample1:group;sample1:group;.."
    groups=[ s.split(":") for s in "{{groups}}".split(";") ]
    groups=pd.DataFrame(groups, columns=["sample","group"] )
    df=df.drop( [ "group" ] , axis=1 )
    df=pd.merge( groups, df, on=["sample"], how="inner" )

df["rep"]=add_rep_suffix( df["group"] )

for experiment in df["Experiment"].tolist() :
    runs=df.loc[ df["Experiment"] == experiment , "Run" ].split(",")
    layout=df.loc[ df["Experiment"] == experiment , "LibraryLayout" ]
    rep=df.loc[ df["Experiment"] == experiment , "rep" ].split(",")
    group=df.loc[ df["Experiment"] == experiment , "group" ]

    if len(runs) == 1 :
        
        run=runs[0]

        if layout == "PAIRED" :
            
            old_name_1=f"{run}_1.fastq.gz" 
            new_name_1=f"{group}.{rep}.read_1.fastq.gz"
            print( f"Renaming {old_name_1} as {new_name_1}" )
            os.rename( old_name_1, new_name_1 )

            old_name_2=f"{run}_2.fastq.gz" 
            new_name_2=f"{group}.{rep}.read_2.fastq.gz"
            print( f"Renaming {old_name_2} as {new_name_2}" )
            os.rename( old_name_2, new_name_2 )

            file = Path(new_name_1)
            new_location = file.parent.parent / file.name  # ../file.txt
            file.rename(new_location)

            file = Path(new_name_2)
            new_location = file.parent.parent / file.name  # ../file.txt
            file.rename(new_location)

        
        
        else:

            old_name=f"{run}.fastq.gz" 
            new_name=f"{group}.{rep}.read_1.fastq.gz"
            print( f"Renaming {old_name} as {new_name}" )
            os.rename( old_name, new_name )

            file = Path(new_name)
            new_location = file.parent.parent / file.name  # ../file.txt
            file.rename(new_location)

    else:    
        
        if layout == "PAIRED" :
            
            new_name_1=f"{group}.{rep}.read_1.fastq.gz"
            new_name_2=f"{group}.{rep}.read_2.fastq.gz"

            files_1=[ f"{s}_1.fastq.gz" for s in runs ]
            files_2=[ f"{s}_2.fastq.gz" for s in runs ]

            files_1_=", ".join(files_1)
            files_2_=", ".join(files_2)

            print( f"Concatenating {files_1_} as {new_name_1}" )
            concat_gz_files( files_1, new_name_1 )
            print( f"Concatenating {files_2_} as {new_name_2}" )
            concat_gz_files( files_2, new_name_2 )

            file = Path(new_name_1)
            new_location = file.parent.parent / file.name  # ../file.txt
            file.rename(new_location)

            file = Path(new_name_2)
            new_location = file.parent.parent / file.name  # ../file.txt
            file.rename(new_location)
        

        else:

            new_name=f"{group}.{rep}.read_1.fastq.gz"
            files=[ f"{s}.fastq.gz" for s in runs ]
            files_=", ".join(files)

            print( f"Concatenating {files_} as {new_name}" )
            concat_gz_files( files, new_name )

            file = Path(new_name)
            new_location = file.parent.parent / file.name  # ../file.txt
            file.rename(new_location)

""",
    desc={
        "downloads_folder": "Downloads folder.",
        "geoacc":"geo accession", 
    },
    container="mpgagebioinformatics/agepy:f92ee75"
)

# sorted.fastq.gz
test_unpigz=jawm.Process( 
    name="test_unpigz",
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], f'{p.var["sraid"]}.fastq'  ) ) ,
    script="""#!/bin/bash
cd {{downloads_folder}}
unpigz {{sraid}}.fastq.gz
""",
    desc={
        "downloads_folder": "Downloads folder.",
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

        input_xlsx=os.path.join( read_geo.var["downloads_folder"], f"{read_geo.var['geoacc']}.samples.xlsx"  )

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

        relabel_geo.execute( )

        jawm.Process.wait( relabel_geo.hash )

    if workflow( "test", workflows ) :

        # for the test workflow we also do something more (just for sra)
        test_unpigz.execute()
        jawm.Process.wait( )
        print("Test completed.")

    sys.exit(0)
