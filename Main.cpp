#include <ctype.h>
#include <time.h>

#include "GenomeSequence.h"
#include "Parameters.h"
#include "StringArray.h"

void printExampleData(FILE *f)
{
	fprintf(f, 
"chr	loc	ref	alleles	snp.Q	av.max.map.Q	depth.cov	NA12891	NA12891.Q	NA12892	NA12892.Q	NA12878	NA12878.Q	hwe	maf	tdt	display\n\
1	52066	t	t/c	64	92.33	36	T/T	29	C/T	65	C/T	64	NO-HWE	MAF-1	1TDT-1	tmp:1:52066\n\
1	91549	a	a/g	76	82.67	89	A/A	79	A/G	71	A/G	67	NO-HWE	MAF-1	1TDT-1	tmp:1:91549\n\
1	223336	c	c/g	100	90.67	70	C/G	100	C/C	56	C/G	100	NO-HWE	MAF-1	1TDT-1	tmp:1:223336\n\
1	695745	g	g/a	97	89	57	A/G	100	G/G	44	A/G	100	NO-HWE	MAF-1	1TDT-1	tmp:1:695745\n\
1	713754	g	g/c	57	89.33	73	C/C	76	C/C	44	C/C	68	NO-HWE	NO-MAF	NOTDT	rs2977670\n\
1	724429	g	g/a	92	90.67	78	G/G	13	A/G	88	A/G	21	LOWQUAL	LOWQUAL	LOWQUAL	tmp:1:724429\n\
1	742429	g	g/a	81	99	61	A/A	86	A/A	46	A/A	85	NO-HWE	NO-MAF	NOTDT	rs3094315\n\
1	742584	a	a/g	100	97.33	86	G/G	82	G/G	83	G/G	100	NO-HWE	NO-MAF	NOTDT	rs3131972\n\
1	744045	a	a/g	71	98	46	G/G	62	G/G	52	G/G	78	NO-HWE	NO-MAF	NOTDT	rs3131969\n\
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
