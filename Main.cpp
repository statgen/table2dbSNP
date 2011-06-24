#include <ctype.h>
#include <time.h>

#include "GenomeSequence.h"
#include "Parameters.h"
#include "StringArray.h"

void printExampleData(FILE *f)
{
    fprintf(f, 
"chr\tloc\tref\talleles\tsnp.Q\tav.max.map.Q\tdepth.cov\tNA12891\tNA12891.Q\tNA12892\tNA12892.Q\tNA12878\tNA12878.Q\thwe\tmaf\ttdt\tdisplay\n\
1\t52066\tt\tt/c\t64\t92.33\t36\tT/T\t29\tC/T\t65\tC/T\t64\tNO-HWE\tMAF-1\t1TDT-1\ttmp:1:52066\n\
1\t91549\ta\ta/g\t76\t82.67\t89\tA/A\t79\tA/G\t71\tA/G\t67\tNO-HWE\tMAF-1\t1TDT-1\ttmp:1:91549\n\
1\t223336\tc\tc/g\t100\t90.67\t70\tC/G\t100\tC/C\t56\tC/G\t100\tNO-HWE\tMAF-1\t1TDT-1\ttmp:1:223336\n\
1\t695745\tg\tg/a\t97\t89\t57\tA/G\t100\tG/G\t44\tA/G\t100\tNO-HWE\tMAF-1\t1TDT-1\ttmp:1:695745\n\
1\t713754\tg\tg/c\t57\t89.33\t73\tC/C\t76\tC/C\t44\tC/C\t68\tNO-HWE\tNO-MAF\tNOTDT\trs2977670\n\
1\t724429\tg\tg/a\t92\t90.67\t78\tG/G\t13\tA/G\t88\tA/G\t21\tLOWQUAL\tLOWQUAL\tLOWQUAL\ttmp:1:724429\n\
1\t742429\tg\tg/a\t81\t99\t61\tA/A\t86\tA/A\t46\tA/A\t85\tNO-HWE\tNO-MAF\tNOTDT\trs3094315\n\
1\t742584\ta\ta/g\t100\t97.33\t86\tG/G\t82\tG/G\t83\tG/G\t100\tNO-HWE\tNO-MAF\tNOTDT\trs3131972\n\
1\t744045\ta\ta/g\t71\t98\t46\tG/G\t62\tG/G\t52\tG/G\t78\tNO-HWE\tNO-MAF\tNOTDT\trs3131969\n\
"
    );
    return;
}

struct arguments {
    arguments() {
        // set argument default values:
        time_t t;
        struct tm *timep;
        char timeBuffer[256];

        t = time(NULL);
        timep = localtime(&t);


        strftime(timeBuffer, sizeof(timeBuffer), "1000G-%Y-%m-%d", timep);
        batch = timeBuffer;

        strftime(timeBuffer, sizeof(timeBuffer), "%Y%m%d%H%M", timep);
        batchID = timeBuffer;

        citation = "SNP's detected by shotgun sequencing a nuclear family with father, mother and child";
        displayFieldName = "display";
        handle = "1000G";
        help = false;
        method = "1000G_TRIOS";
        moltype = "Genomic";
        noHeader = false;
        organism = "Homo sapiens";
        printExampleData = false;
        sampleSize = 4;
        theReference = "/home/1000G/data/chromosomes5000.fa";
    }
    void writeHeader(FILE *);
    void usage(int argc, char **argv);
    String  batch;      // used in header
    String  batchID;    // used for ACCESSION value in per SNP report
    String  citation;
    String  comment;
    String  displayFieldName;
    String  handle;
    bool    help;
    String  method;
    String  moltype;
    bool    noHeader;
    String  organism;
    bool    printExampleData;
    int     sampleSize;
    String  theReference;
};

void arguments::writeHeader(FILE *f)
{
    if(noHeader) return;
    fprintf(f, "TYPE: SNPASSAY\n");
    fprintf(f, "HANDLE: %s\n", (const char *) handle);
    fprintf(f, "BATCH: %s\n", (const char *) batch);
    fprintf(f, "MOLTYPE: %s\n", (const char *) moltype);
    fprintf(f, "METHOD: %s\n", (const char *) method);
    fprintf(f, "SAMPLESIZE: %d\n", sampleSize);
    fprintf(f, "ORGANISM:%s\n", (const char *) organism);
    fprintf(f, "CITATION:%s\n", (const char *) citation);
    if(comment!="") fprintf(f, "COMMENT:%s\n", (const char *) comment);
    fprintf(f, "||\n");
}

void arguments::usage(int argc, char **argv)
{   
    fprintf(stderr,"usage: %s [options]\n", argv[0]);
    fprintf(stderr, "\t--batch\n");
    fprintf(stderr, "\t--batchID\n");
    fprintf(stderr, "\t--comment\n");
    fprintf(stderr, "\t--handle\n");
    fprintf(stderr, "\t--help\n");
    fprintf(stderr, "\t--moltype\n");
    fprintf(stderr, "\t--noHeader\n");
    fprintf(stderr, "\t--printExampleData -> prints out an example of the expected data format\n");
    fprintf(stderr, "\t--sampleSize\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s reads stdin, and writes to stdout.\n", argv[0]);
    fprintf(stderr, "\n");
    fprintf(stderr, "example:\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "%s --printExampleData | %s >result.txt\n", argv[0], argv[0]);
    exit(1);
}

int main (int argc, char ** argv)
{
    struct arguments args;
    ParameterList pl;
    int argsUsed;

    int snpLinkFieldIndex;

    StringArray columnHeaders;
    StringArray inputFields;
    String inputLine;

    BEGIN_LONG_PARAMETERS(parameters)
    LONG_PARAMETER_GROUP("Input Files")
        LONG_STRINGPARAMETER("reference", &args.theReference)
    LONG_PARAMETER_GROUP("Options")
        LONG_STRINGPARAMETER("batch", &args.batch)
        LONG_STRINGPARAMETER("batchID", &args.batchID)
        LONG_STRINGPARAMETER("citation", &args.citation)
        LONG_STRINGPARAMETER("comment", &args.comment)
        LONG_STRINGPARAMETER("displayFieldName", &args.displayFieldName)
        LONG_PARAMETER("printExampleData", &args.printExampleData)
        LONG_STRINGPARAMETER("handle", &args.handle)
        LONG_PARAMETER("help", &args.help)
        LONG_STRINGPARAMETER("moltype", &args.moltype)
        LONG_PARAMETER("noHeader", &args.noHeader)
        LONG_INTPARAMETER("sampleSize", &args.sampleSize)
    END_LONG_PARAMETERS();

    pl.Add(new LongParameters("Available Options", parameters));
    argsUsed = pl.ReadWithTrailer(argc, argv);

    if(args.help) {
        pl.Status();
        args.usage(argc, argv);
    }

    if(args.printExampleData) {
        printExampleData(stdout);
        exit(0);
    }


    GenomeSequence reference;

    if(reference.setReferenceName(args.theReference.c_str())) {
        fprintf(stderr,"failed to open reference file %s.\n",
        (const char *) args.theReference);
        exit(1);
    }

    if(reference.open()) {
        exit(1);
    }

    args.writeHeader(stdout);
    int lineNumber = 0;

    while(!feof(stdin)) {
        genomeIndex_t genomeIndex;
        int chromosome;

        inputLine.ReadLine(stdin);
        inputFields.ReplaceColumns(inputLine);
        if(lineNumber++ == 0) {
            columnHeaders = inputFields;
            snpLinkFieldIndex = -1;
            for(int i = 0; i<inputFields.Length(); i++) {
                if(inputFields[i]==args.displayFieldName) {
                    snpLinkFieldIndex = i;
                    break;
                }
            }
            if(snpLinkFieldIndex < 0) {
                fprintf(stderr,"unable to find column header named '%s'!\n", args.displayFieldName.c_str());
                exit(1);
            }
            continue;
        }
#define REQUIRED_ARG_COUNT 4
        if(inputFields.Length() < REQUIRED_ARG_COUNT) {
            fprintf(stderr,"%d: mal formed input line - must have at least %d fields\n",
                lineNumber,
                REQUIRED_ARG_COUNT);
            fprintf(stderr,"%d: %s\n", lineNumber, (const char *) inputLine);
            continue;
        }
        String chromosomeName = inputFields[0];
        unsigned int chromosomeIndex = atoi(inputFields[1]);
        String original = inputFields[2];
        String substitute = inputFields[3];

#if defined(DEBUG)
        printf("---------------------------\n");
        printf("chromosome='%s', chromosomeIndex='%u', original='%s', sub='%s'\n",
            (const char *) chromosomeName, chromosomeIndex, (const char *) original, (const char *) substitute);

#endif

        chromosome = reference.getChromosome(chromosomeName.c_str());
        genomeIndex = reference.getGenomePosition(chromosome, chromosomeIndex);
        
        //
        // Sanity check: the genome Index should be valid:
        //
        if(genomeIndex == INVALID_GENOME_INDEX) {
            fprintf(stderr, "%d: chromosome %s, position %u was not found in the reference genome.\n",
                lineNumber,
                (const char *) chromosomeName,
                chromosomeIndex);
            continue;
        }

        //
        // Sanity check: the genome Index should be correctly mapped:
        //
        if(genomeIndex > reference.sequenceLength()) {
            fprintf(stderr, "%d: chromosome %s, position %u exceeds length of reference genome - skipped.\n",
                lineNumber,
                (const char *) chromosomeName,
                chromosomeIndex);
            continue;
        }

        //
        // Sanity check: the reference Genome and the claimed "original" should be the same
        // base.
        //
        if(tolower(reference[genomeIndex]) != tolower(original[0])) {
            fprintf(stderr, "%d: chromosome %s, position %u, expected %c, but have %c - skipped.\n",
                lineNumber,
                chromosomeName.c_str(),
                chromosomeIndex,
                original[0],
                reference[genomeIndex]);
            continue;
        }
#define LEFT_FLANK_LEN 200
#define RIGHT_FLANK_LEN 200

        // Assumes that lineNumber == 2 is the start of the data.
        // We use underscores at the request of Lon Phan, who finds
        // them easier to parse for their purposes.
        printf("SNP:%s|%s_%d_chr%s_%u\n", args.handle.c_str(), args.batchID.c_str(), lineNumber - 1, chromosomeName.c_str(), chromosomeIndex);
        if(strncmp(inputFields[snpLinkFieldIndex].c_str(), "rs", 2)==0) printf("SNP_LINK:NCBI|%s\n",
                inputFields[snpLinkFieldIndex].c_str());
        printf("ACCESSION:NC_%06d\n", chromosome+1);  // 
        printf("LOCATION:%u\n", chromosomeIndex);       // XXX relative to ACCESSION?
        printf("SAMPLESIZE:%d\n",args.sampleSize);
        printf("LENGTH:%d\n",LEFT_FLANK_LEN + RIGHT_FLANK_LEN + 1);

        genomeIndex_t i;

        printf("5'_FLANK:");

        for(i=genomeIndex-200;i<genomeIndex;i++) printf("%c", reference[i]);
        printf("\n");

        printf("OBSERVED:%s\n", substitute.ToUpper().c_str());

        printf("3'_FLANK:");
        for(i=genomeIndex+1;i<genomeIndex+201;i++) printf("%c", reference[i]);
        printf("\n");

        printf("COMMENT:\n");
        for(int i=REQUIRED_ARG_COUNT; i<inputFields.Length(); i++) {
            printf("%s:%s\n", (const char *) columnHeaders[i], (const char *) inputFields[i]);
        }
        printf("||\n");

    }
}
