#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <ctype.h>
#define MAXSTR 51200
int count_headers(FILE *);

/* revert_headers: FastaAlternateReferenceMaker changes the chromosome names; this changes
them back to match the reference genome */
/*********************************************************************************************/
void main(int argc,char **argv)
{
        FILE *f, *g;
	char **header = NULL, line[MAXSTR];
	char output_filename[MAXSTR];
	int header_cnt = 0, cnt = 0, i = 0;
	
	if (argc < 3) {
		printf("Usage: %s <reference_genome> <genome_to_fix> [output_prefix]\n", argv[0]);
		printf("  reference_genome: Original reference genome file\n");
		printf("  genome_to_fix:    Genome file with modified headers to fix\n");
		printf("  output_prefix:    Optional output filename (default: tmpc.fasta)\n");
		exit(0);
	}
	
	// Set output filename - use provided prefix or default
	if (argc >= 4) {
		// Check if the provided name already has .fasta extension
		if (strstr(argv[3], ".fasta") != NULL || strstr(argv[3], ".fa") != NULL) {
			strcpy(output_filename, argv[3]);
		} else {
			snprintf(output_filename, MAXSTR-1, "%s.fasta", argv[3]);
		}
	} else {
		strcpy(output_filename, "tmpc.fasta");
	}
	
	printf("Reverting headers from %s to headers from %s\n", argv[2], argv[1]);
	printf("Output will be written to: %s\n", output_filename);
	
	if ((f = fopen(argv[1],"r")) == NULL) {
                printf("Can't open reference file %s\n",argv[1]);
                exit(1);
        }
	if ((g = fopen(output_filename,"w")) == NULL) {
                printf("Can't open output file %s\n", output_filename);
                fclose(f);
                exit(1);
        }
	
	header_cnt = count_headers(f);
	printf("Found %d headers in reference genome\n", header_cnt);
	
	if (header_cnt == 0) {
		printf("No headers found in reference genome\n");
		fclose(f);
		fclose(g);
		exit(1);
	}
	
	if ((header=calloc(header_cnt, sizeof(char *)))==NULL) {
                printf("Memory allocation error for headers\n");
                fclose(f);
                fclose(g);
                exit(1);
        }
        
	// Read headers from reference genome
	i = 0;
	while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>' && i < header_cnt) {
			if ((header[i]=calloc((strlen(line) + 1), sizeof(char)))==NULL) {
                                printf("Memory allocation error for header %d\n", i);
                                exit(1);
                        }
                        strcpy(header[i], line);
                        i++;
		}
	}
	fclose(f);
	
	// Process the genome file to fix
	if ((f = fopen(argv[2],"r")) == NULL) {
                printf("Can't open genome file to fix: %s\n",argv[2]);
                exit(1);
        }
        
	i = 0;
	while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>') {
			cnt++;
			if (cnt > header_cnt) {
				printf("ERROR: Genome file has more headers (%d) than reference (%d)\n", 
				       cnt, header_cnt);
				exit(1);
			}
			fprintf(g, "%s", header[i]);
			printf("Replaced: %s", line);
			printf("With:     %s", header[i]);
			i++;
		}
		else {
			fprintf(g, "%s", line);
		}
	}
	
	if (cnt < header_cnt) {
		printf("WARNING: Genome file has fewer headers (%d) than reference (%d)\n", 
		       cnt, header_cnt);
	}
	
	printf("Successfully processed %d sequences\n", cnt);
	printf("Output written to: %s\n", output_filename);
	
	// Cleanup
	for (i = 0; i < header_cnt; i++) {
		if (header[i] != NULL) 
			free(header[i]);
	}
	if (header != NULL) 
		free(header);
		
	fclose(f);
	fclose(g);
	exit(0);
}

/*********************************************************************************************/
int count_headers(FILE *f)
{
        int cnt = 0;
        char line[MAXSTR];
        while(fgets(line,MAXSTR - 1,f) != NULL) {
                if (line[0] == '>')
                        cnt++;
        }
        rewind(f);
        return(cnt);
}
/*********************************************************************************************/