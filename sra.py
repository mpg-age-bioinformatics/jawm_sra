import jawm
import os

prefetch=jawm.Process( 
    name="prefetch",
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], f'{p.var["sraid"]}.fastq.gz'  ) ) ,
    script="""#!/bin/bash
cd {{downloads_folder}}
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
    when=lambda p: not os.path.isfile( os.path.join( p.var["downloads_folder"], f'{p.var["sraid"]}.fastq.gz'  ) ) ,
    script="""#!/bin/bash
cd {{downloads_folder}}

if [ -f "{{sraid}}.sralite" ]; then
    fastq-dump --split-3 "./{{sraid}}.sralite"
    rm -rf "{{sraid}}.sralite"
else
    fastq-dump --split-3 "./{{sraid}}"
    rm -rf "{{sraid}}"
fi

if [ -f {{sraid}}.fastq ] ; then pigz {{sraid}}.fastq ; fi
if [ -f {{sraid}}_1.fastq ] ; then pigz {{sraid}}_1.fastq ; fi
if [ -f {{sraid}}_2.fastq ] ; then pigz {{sraid}}_2.fastq ; fi

rm -rf {{sraid}}
""",
    desc={
        "downloads_folder": "Downloads folder.",
        "sraid":"sra_id", 
    },
    container="mpgagebioinformatics/sra:3.2.1"
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

    workflows, var, args, unknown_args=jawm.utils.parse_arguments( ["main","sra","test"] )

    # usage: 
    if workflow( ["main","sra","test"], workflows ) :

        # pre-download images to avoid parallel 
        # processes concurring to download the same image
        images=get_image( prefetch.container )

        sra_jobs=[]
        for sraid in prefetch.var["sraid"].split(",") : 
            
            prefetch_=prefetch.clone()
            prefetch_.var["sraid"]=sraid
            prefetch_.execute()

            fastq_dump_=fastq_dump.clone()
            fastq_dump_.execute( prefetch_.hash )

            sra_jobs.append( fastq_dump_.hash )
        
        jawm.Process.wait( sra_jobs )

    if workflow( "test", workflows ) :

        # for the test workflow we also do something more (just for sra)
        test_unpigz.execute()
        jawm.Process.wait( )
        print("Test completed.")

    sys.exit(0)
