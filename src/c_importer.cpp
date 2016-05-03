#include <Rcpp.h>
#include <fstream>
#include <math.h>
#include <algorithm>
#include <vector> 

using namespace Rcpp;
using namespace std;
#define HEADERSIZE 6
#define INS_IONTOF "iontof"
#define CEIL_MZ 5000

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
List read_GRD(char* filename, int imageSize);
List read_ITZIP(char* scans, 
                char* tofs, 
                char* shots, char* coords, int imageSize);
List read_BIF(std::string file, std::string instrument);
List calibrate_ITZIP(List uncalibrate, float upperMass);
IntegerVector verify_input_ITZIP(CharacterVector inputs);
int is_existed(CharacterVector inputs, char *fileType);
float round_with_precision(float value, int precision);

struct c_unique {
        int current;
        c_unique() {current=0;}
        int operator()() {return ++current;}
} UniqueNumber;

//' @name import
//' @description import is the C++ code for importing iontof raw data
//' @title import is the C++ code for importing iontof raw data
//' @param rFilename CharacterVector
//' @param fType CharacterVector
//' @param imageSize int
//' @param upperMass int
//' @return imported binary raw data
// [[Rcpp::export]]
List import(CharacterVector rFilename, 
            CharacterVector fType, int imageSize, float upperMass){
        //printf("\nImport ");
        char* filename = rFilename[0];
        char* filetype = fType[0];
        List empty(0);
        if(strcmp(filetype, "grd") == 0){
                //printf("GRD format.\n");
                return(calibrate_ITZIP(
                        read_GRD(filename, imageSize), 
                        upperMass));
        } else if(strcmp(filetype, "itzip") == 0){
                //printf("ITZIP format.\n");
                IntegerVector order = verify_input_ITZIP(rFilename);
        if(sum(order) >= 6){
            char *scans = rFilename[order[0]];
            char *tofs = rFilename[order[1]];
            char *shots = rFilename[order[2]];
            char *coords = rFilename[order[3]];
            return(calibrate_ITZIP(
                    read_ITZIP(scans, tofs, shots, coords, imageSize), 
                    upperMass));
//            return(read_ITZIP(scans, tofs, shots, coords, imageSize));
        } else {
//printf("unknown format.\n");
//printf("At the moment we just support \"itzip\" and \"grd\" format.\n");
        }
    }
    return(empty);
}

// \\[[Rcpp::export]]
//NumericMatrix build_mass_spectra(NumericMatrix importer_out, int imageSize){
//    int tof_max = importer_out(0,1);
//    int tof_min = importer_out(0,1);
//    int count = 1;
//    int count_ref = 2;
//    int total_ion_count = importer_out.nrow();
//    NumericVector tofs = unique(importer_out(_,1));
//    int pixels = imageSize * imageSize;
//    NumericMatrix mass_spectra(pixels, tofs.size());
//    
//    return(importer_out);
//}

List read_GRD(char* filename, int imageSize){
    FILE *f = fopen(filename, "rb");
    if(f == NULL){
        perror("Error");
        List dataValues(0);
        return(dataValues);
    }

	fseek(f, 0, SEEK_END);
	long long fileLen = ftell(f);
	fseek(f, 0, SEEK_SET);
    long long dataBlock = fileLen / 20;
//    printf("Ion counts : %d\n", dataBlock);
    //cout << "Ion counts : " << dataBlock << endl;
    
    NumericMatrix dataValues(dataBlock, 4);
    long count = 0;
    
    unsigned int dataValue[5];
//    imageSize = 1;
    size_t read_res = 1;
    //printf("start to read data");
    while(!feof(f)){
        read_res = read_res && fread(dataValue, sizeof(unsigned int), 5, f);
        
        if(read_res == 1){
            // x * y
            dataValues(count, 0) = dataValue[2] * imageSize + dataValue[3] + 1;
            // tofs
            dataValues(count, 1) = dataValue[4];
            // scans
            dataValues(count, 2) = dataValue[0];
            // shots
            dataValues(count, 3) = dataValue[1];
            count++;
        }
    }
    //printf("\nDone.");
    fclose(f);
    List data(0);
    int highestTofs = dataValues(which_max(dataValues(_,1)), 1);
    data["highestValues"] = highestTofs;
    data["importedMatrix"] = dataValues;
    return(data);
}

// ITZIP FORMAT
// The inputs must contains 4 files : scans, tofs, shots and coords.
// Output : the order of : scans, tofs, shots, coords.
IntegerVector verify_input_ITZIP(CharacterVector inputs){
    //printf("Verify inputs.\n");
    IntegerVector order(4);
    char scan_ext[] = ".scans";
    char tof_ext[] = ".tofs";
    char shot_ext[] = ".shots";
    char coord_ext[] = ".coords";
    order[0] = is_existed(inputs, scan_ext);
    order[1] = is_existed(inputs, tof_ext);
    order[2] = is_existed(inputs, shot_ext);
    order[3] = is_existed(inputs, coord_ext);
    
    if(sum(order) >= 6){
        //printf("Verify passed.\n");
    } else {
        //printf("Verify failed.\n");
    }
    return(order);
}

// Check whether the required file is existed.
int is_existed(CharacterVector inputs, char *fileType){
    for(int i = 0; i < inputs.size(); i++){
        if(strstr(inputs[i], fileType) != NULL){
            //printf(" - The inputs contains %s file.\n", fileType);
            return(i);
        }
    }
    //printf(" - The inputs lack of %s file.\n", fileType);
    return(-1);
}

float round_with_precision(float value, int precision){
    long precision_base = pow(10, precision);
    return roundf(value * precision_base) / precision_base;
}

// This function is used for read and combine scans, tofs, shots and coords 
// data.
List read_ITZIP(char* scans, 
                char* tofs, 
                char* shots, 
                char* coords, 
                int imageSize){
    FILE *fscans = fopen(scans, "rb");
    FILE *ftofs = fopen(tofs, "rb");
    FILE *fshots = fopen(shots, "rb");
    FILE *fcoords = fopen(coords, "rb");
    List data(0);
    
    if(fscans == NULL || ftofs == NULL || fshots == NULL || fcoords == NULL){
        perror("\nError");
        NumericMatrix dataValues(0,0);
        data["importedMatrix"] = dataValues;
        return(data);
    }

    fseek(fcoords, 0, SEEK_END);
    // Divide by 4 because they used 4 bytes for each data point
    // Divide by 2 because each coord contains 2 values (x and y)
	int dataBlock = ftell(fcoords) / 4 / 2;
	fseek(fcoords, 0, SEEK_SET);
    //printf("\nIon counts : %d", dataBlock);
    
    NumericMatrix dataValues(dataBlock, 4);
    int count = 0;
    
    unsigned int scan_data[1], tof_data[1], shot_data[1], coord_data[2];
//    imageSize -= 1;
    size_t read_res = 1;
    while(!feof(fscans)){
        read_res = read_res && 
            fread(scan_data, sizeof(unsigned int), 1, fscans);
        read_res = read_res && 
            fread(tof_data, sizeof(unsigned int), 1, ftofs);
        read_res = read_res && 
            fread(shot_data, sizeof(unsigned int), 1, fshots);
        read_res = read_res && 
            fread(coord_data, sizeof(unsigned int), 2, fcoords);
        if(read_res == 1){
            // x * y
            dataValues(count, 0) = coord_data[0] * imageSize + 
                coord_data[1] + 1;
            // TOF data
            dataValues(count, 1) = tof_data[0];
            // Scan data
            dataValues(count, 2) = scan_data[0];
            // Shooting time
            dataValues(count, 3) = shot_data[0];
            count++;
        }
    }
    fclose(fscans);
    fclose(ftofs);
    fclose(fshots);
    fclose(fcoords);
    //printf("\nDone.");
    int highestTofs = dataValues(which_max(dataValues(_,1)), 1);
    data["highestTofs"] = highestTofs;
    data["importedMatrix"] = dataValues;
    return(data);
}

//This function is used to calibrate ITZIP format.
// [[Rcpp::export]]
List calibrate_ITZIP(List uncalibrate, float upperMass){
    //printf("\nCalibrate data.");
    NumericMatrix uncalMatrix = uncalibrate["importedMatrix"];
    int highestTofs = uncalMatrix(which_max(uncalMatrix(_,1)), 1);
    uncalMatrix(_,1) = pow((uncalMatrix(_,1) / 
                                (highestTofs / sqrt(upperMass))), 2);
    //printf("\nDone calibrating.\n");
    List calibrated(0);
    calibrated["highestTofs"] = highestTofs;
    calibrated["calibratedMatrix"] = uncalMatrix;
    return(calibrated);
}

// [[Rcpp::export]]
List read_BIF(std::string file, std::string instrument) { 
    List bif;
    int found = file.find_first_of("~");
    if(found != -1)
        file = file.replace(0,1,getenv("HOME"));
    
    ifstream fbin (file.c_str(), ios::binary | ios::in);
    if (!fbin) {
        //cout << "Could not open file " << file;
        return bif;
    }
    
    // Skip the header for iontof format
    if(instrument.compare(INS_IONTOF) == 0){
        fbin.seekg (0, ios::beg);
        fbin.seekg(HEADERSIZE * sizeof(char));
    }
    
    uint16_t nIntervals, nXPixels, nYPixels, nStepSize;
    int nPixels, nData;
    
    fbin.read((char *)&nIntervals, sizeof(nIntervals));
    fbin.read((char *)&nXPixels, sizeof(nXPixels));
    fbin.read((char *)&nYPixels, sizeof(nYPixels));
    nPixels = nXPixels * nYPixels;
    nData = nPixels * nIntervals;
    
    vector<int> dim1Flip(nXPixels), dim2Flip(nYPixels);
    vector<int> ids;
    vector<uint16_t> imageVec;
    vector<float> middles, lowers, uppers, massValues;
//    vector<>
    
    std::generate(dim1Flip.begin(), dim1Flip.end(), UniqueNumber);
    std::generate(dim2Flip.begin(), dim2Flip.end(), UniqueNumber);
    
    if(instrument.compare(INS_IONTOF) == 0){
        nStepSize = 1;
        // reverese dim2Flip for IONTOF
        std::reverse(dim2Flip.begin(),dim2Flip.end());
        
    } else {
        nStepSize = 2;
        // reverse dim1Flip and dim2Flip for UlvacPhi
        std::reverse(dim2Flip.begin(),dim2Flip.end());
        std::reverse(dim1Flip.begin(),dim1Flip.end());
    
    }
    
    //cout << "importing" << " " << file << endl;
    //cout << nIntervals << " m/z\'s" << endl;
    //cout << nXPixels << " x pixels" << endl;
    //cout << nYPixels << " y pixels" << endl;
    //cout << nXPixels * nYPixels << " pixels" << endl;
    //cout << nData << " data values" << endl;
    //cout << "reading binary data..." << endl;
    
    uint32_t id, imageDatum;
    float lower, middle, upper;
    int count = 0;
    
    for(int i = 0; i < nIntervals; i++){
        fbin.read((char *)&id, sizeof(uint32_t) / nStepSize);
        fbin.read((char *)&lower, sizeof(float));
        fbin.read((char *)&middle, sizeof(float));
        fbin.read((char *)&upper, sizeof(float));
        
        if(lower > CEIL_MZ) { 
            count++;
            for(int j = 0; j < nPixels; j++)
                fbin.read((char *)&imageDatum, sizeof(uint32_t) / nStepSize);
            continue;
        }
        
        for(int j = 0; j < nPixels; j++){
            fbin.read((char *)&imageDatum, sizeof(uint32_t) / nStepSize);
            imageVec.push_back(imageDatum);
        }
        
        ids.push_back(id);
        lowers.push_back(round_with_precision(lower, 6));
        middles.push_back(round_with_precision(middle, 6));
        massValues.push_back(round_with_precision(middle, 6));
        uppers.push_back(round_with_precision(upper, 6));
    }
    fbin.close();
    nIntervals = nIntervals - count;
    nData = nIntervals * nPixels;
    IntegerVector imageData(nData);
    IntegerVector bifDim = IntegerVector::create(nXPixels, 
                                                 nYPixels, 
                                                 nIntervals);
    for(int i = 0; i < nData; i++){
        imageData(i) = imageVec[i];
    }
    //cout << "reshaping image data to 3D array..." << endl;
    imageData.attr("dim") = bifDim;
    
    bif["ids"] = ids;
    bif["lowers"] = lowers;
    bif["middles"] = middles;
    bif["uppers"] = uppers;
    bif["nPixels"] = nPixels;
    bif["nXPixels"] = nXPixels;
    bif["nYPixels"] = nYPixels;
    bif["nIntervals"] = nIntervals;
    bif["nData"] = nData;
    bif["dim1Flip"] = dim1Flip;
    bif["dim2Flip"] = dim2Flip;
    bif["massValues"] = massValues;
    bif["imageData"] = imageData;
    return bif;
}




