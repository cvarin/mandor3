//contains utilities for Mandor postprocessor such as folder cleaning and so forth

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<stdarg.h>
#include"utils.h"
#include"params.h"
#include"types.h"

void cfg_read_str(char *keyword, char* value, char *cfgfile){
	//looks for the keyword in the cfg file and returns its string value

	int j, nfounds, len;
	char str[MAXSTRLEN];			//temporary storage for the string read
	
	nfounds = 0;				//number of the keywords found in the cfg file
	len = strlen(keyword);			//length of the searched string
	
	FILE *fp;
	OPEN_FILE(fp, cfgfile, "rt");
	
	for(j = 0; j < MAXSTRLEN; j++) str[j] = '\0';
	//scanning the file
	while(fgets(str, MAXSTRLEN, fp)){
		if(str[0] != '#'){
			int i = 0;
			while(str[i] == keyword[i]) i++;
			if((str[i] == ' ')&&(str[++i] == '=')&&(str[++i] == ' ')&&(i == (len+2))){
				nfounds++;
				int k = 0;
				while((str[++i] != '\n')&&(str[i] != '\0')) value[k++] = str[i];
				value[k] = '\0';
			}
		}
	}
	if(nfounds == 0) error(FUNCTION, "the option %s cannot be found in the cfg file", keyword);
	if(nfounds > 1) error(FUNCTION, "multiple definition of the option %s in the cfg file", keyword);
	fclose(fp);
}

int cfg_read_int(char *keyword, char *cfgfile){
	//looks for the keyword in the cfg file and returns its int value

	int j, nfounds, len, val;
	char str[MAXSTRLEN];			//temporary storage for the string read
	char tmp[STRLEN];
	
	nfounds = 0;				//number of the keywords found in the cfg file
	len = strlen(keyword);			//length of the searched string
	
	FILE *fp;
	OPEN_FILE(fp, cfgfile, "rt");
	
	for(j = 0; j < MAXSTRLEN; j++) str[j] = '\0';
	//scanning the file
	while(fgets(str, MAXSTRLEN, fp)){
		if(str[0] != '#'){
			int i = 0;
			while(str[i] == keyword[i]) i++;
			if((str[i] == ' ')&&(str[++i] == '=')&&(str[++i] == ' ')&&(i == (len+2))){
				nfounds++;
				sscanf(str, "%s = %d", tmp, &val);
			}
		}
	}
	if(nfounds == 0) error(FUNCTION, "the option %s cannot be found in the cfg file", keyword);
	if(nfounds > 1) error(FUNCTION, "multiple definition of the option %s in the cfg file", keyword);
	fclose(fp);
	return(val);
}

double cfg_read_double(char *keyword, char *cfgfile){
	//looks for the keyword in the cfg file and returns its double value

	int j, nfounds, len;
	float val;
	char str[MAXSTRLEN];			//temporary storage for the string read
	char tmp[STRLEN];
	
	nfounds = 0;				//number of the keywords found in the cfg file
	len = strlen(keyword);			//length of the searched string
	
	FILE *fp;
	OPEN_FILE(fp, cfgfile, "rt");
	
	for(j = 0; j < MAXSTRLEN; j++) str[j] = '\0';
	//scanning the file
	while(fgets(str, MAXSTRLEN, fp)){
		if(str[0] != '#'){
			int i = 0;
			while(str[i] == keyword[i]) i++;
			if((str[i] == ' ')&&(str[++i] == '=')&&(str[++i] == ' ')&&(i == (len+2))){
				nfounds++;
				sscanf(str, "%s = %f", tmp, &val);
			}
		}
	}
	if(nfounds == 0) error(FUNCTION, "the option %s cannot be found in the cfg file", keyword);
	if(nfounds > 1) error(FUNCTION, "multiple definition of the option %s in the cfg file", keyword);
	fclose(fp);
	return((double)val);
}

void preparefolders(){
	//empties and/or (if necessary) creates folders dat and pics
	system("rm -rf pics/*");
	system("rm -rf dat/*");
	system("mkdir -p pics");
	system("mkdir -p dat");
}











