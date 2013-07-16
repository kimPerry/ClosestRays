//============================================================================
// Name        : graphing.cpp
// Author      : Kim Perry
// Description : Plotting points without POCAs
//============================================================================

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <cmath>

using namespace std;

struct muonTrack
{
	double time;
	double inDirX;
	double inDirY;
	double inDirZ;
	double inPointX;
	double inPointY;
	double inPointZ;
	double outDirX;
	double outDirY;
	double outDirZ;
	double outPointX;
	double outPointY;
	double outPointZ;
	double momentumIn;
	double momentumOut;
	double scatAngleX;
	double scatAngleY;
	double scatAngleTotal;
	double pocaX;
	double pocaY;
	double pocaZ;
};

struct pocas
{
	vector<double> pocaX;
	vector<double> pocaY;
	vector<double> pocaZ;
	vector<double> scatteringAngle;
};

int fileReadIn(ifstream & fin, vector<muonTrack>& muonTracks);
void checkForWhiteSpace(vector<double>lowX, vector<double>lowY, vector<double>lowZ, vector<double>& highX, vector<double>& highY, vector<double>& highZ, vector<double>& scat);
void initializeLookups(int * lowScatlookupX, int * lowScatlookupY, int * lowScatlookupZ);

int main() {

	char filename[] = "MTTrackEpoch_852.csv";
    int unScatPointsPlotted = 0;
	int highScatterCountBefore = 0;

	vector<double> lowScatlookupX;
	vector<double> lowScatlookupY;
	vector<double> lowScatlookupZ;

	ifstream fin;
	ofstream fout;

	//initializeLookups(lowScatlookupX, lowScatlookupY, lowScatlookupZ, lookupLength);

	// array of muon track data
	vector<muonTrack> muonTracks;

	vector<int>::iterator it;

	pocas * highScatterPocas = new pocas;

	// read in file
	cout << "Reading file..." << endl;
	fin.open(filename);
	fileReadIn(fin, muonTracks);
	fin.close();

        int lineCount = muonTracks.size();
        cout << lineCount << " tracks" << endl;
		cout << muonTracks[lineCount-1].pocaX << endl;
		cout << muonTracks[lineCount-2].pocaX << endl;

	// graph stuff
	for(int i = 0; i < lineCount; i++)
	{
		// if track has valid POCA values
		if(muonTracks[i].pocaX > 0 && muonTracks[i].pocaY > 0 && muonTracks[i].pocaZ > 0)
		{
			// if track has scattering angle over a certain threshold
			if(muonTracks[i].scatAngleTotal > .35)
			{
				highScatterPocas->pocaX.push_back(muonTracks[i].pocaX);
				highScatterPocas->pocaY.push_back(muonTracks[i].pocaY);
				highScatterPocas->pocaZ.push_back(muonTracks[i].pocaZ);
				highScatterPocas->scatteringAngle.push_back(muonTracks[i].scatAngleTotal);
				highScatterCountBefore++;
			}

			// else add to low scatter list
			else
			{
				if(muonTracks[i].scatAngleTotal < .002)
				{
					lowScatlookupX.push_back(muonTracks[i].pocaX);
					lowScatlookupY.push_back(muonTracks[i].pocaY);
					lowScatlookupZ.push_back(muonTracks[i].pocaZ);
					unScatPointsPlotted++;
				}
			}
		}
	}


    cout << "high scattering points before checking: " << highScatterCountBefore << endl;
    cout << "low scattering points: " << unScatPointsPlotted << endl;
	cout << "Looking for white space..." << endl;

	checkForWhiteSpace(lowScatlookupX, lowScatlookupY, lowScatlookupZ, highScatterPocas->pocaX, highScatterPocas->pocaY,
			           highScatterPocas->pocaZ, highScatterPocas->scatteringAngle);


	// read out to file
	cout << "Writing file..." << endl;
	fout.open("data.dat");
	cout << "file: " << int(highScatterPocas->pocaX.size()) << endl;
	for(int i = 0; i < int(highScatterPocas->pocaX.size()); i++)
	{
		fout << highScatterPocas->pocaX[i] << "\t" << highScatterPocas->pocaY[i] << "\t" << highScatterPocas->pocaZ[i] << "\t" << highScatterPocas->scatteringAngle[i] * 100 << endl;
	}
	fout.close();

	cout << "high scattering points after checking: " << int(highScatterPocas->pocaX.size()) << endl;

	return 0;
}

void checkForWhiteSpace(vector<double>lowX, vector<double>lowY, vector<double>lowZ, vector<double>& highX, vector<double>& highY, vector<double>& highZ,
						vector<double>& scat)
{

	vector<int> yIndices;
	vector<int> zIndices;
	vector<int> badIndices;

	double percentage = .05;

	// loop through all x's
	//cout << int(lowX.size()) << endl;
	for(int i = 0; i < int(highX.size()); i++)
	{
		for( int j = 0; j < int(lowX.size()); j++)
		{
			if((lowX[j] - (lowX[j] * percentage)) < highX[i] && highX[i] < (lowX[j] + (lowX[j] * percentage)))
			{
				// for every matching X add that index to a list to check the y's of
				yIndices.push_back(j);
				//cout << "low: " << lowX[j] << " high: " << highX[i] << endl;
			}
		}

		//cout << int(yIndices.size()) << endl;
		for( int j = 0; j < int(yIndices.size()); j++)
		{

			if((lowY[yIndices[j]] - (lowY[yIndices[j]] * percentage)) < highY[i] && highY[i] < (lowY[yIndices[j]] + (lowY[yIndices[j]] * percentage)))
			{
				//cout << "low: " << lowY[yIndices[j]] << " high: " << highY[i] << endl;
				// for every matching Y add that index to a list to check the z's of
				zIndices.push_back(j);
			}
		}

		for(int j = 0; j < int(zIndices.size()); j++)
		{
			// see if high Z POCA is within +- 10% of a low Z POCA
			if((lowZ[zIndices[j]] - (lowZ[zIndices[j]] * percentage)) < highZ[i] && highZ[i] < (lowZ[zIndices[j]] + (lowZ[zIndices[j]] * percentage)))
			{
				//cout << "low: " << lowZ[zIndices[j]] << " high: " << highZ[i] << endl;
				// if found remove from high scatter list
				badIndices.push_back(i);
			}
		}
		//cout << "x's: " << lowX.size() <<  " y's: " << yIndices.size() << " z's: " << zIndices.size() << " bad: " << badIndices.size() << endl;
		yIndices.clear();
		zIndices.clear();
	}

	cout << "bad indices found: " << int(badIndices.size()) << endl;

	for(int i = int(badIndices.size()) -1; i > -1; i--)
	{
		//cout << badIndices[i] << endl;

		highX.erase(highX.begin() + badIndices[i]);
		highY.erase(highY.begin() + badIndices[i]);
	    highZ.erase(highZ.begin() + badIndices[i]);
		scat.erase(scat.begin() + badIndices[i]);
	}
cout << "new size: " << highX.size() << endl;
cout << "finished function" << endl;
	/*
		index = it - lowX.begin();
*/
}

int fileReadIn(ifstream & fin, vector<muonTrack> & muonTracks)
{
	string temp;

    muonTrack dummyMuon;

	int index = 0;

	muonTracks.push_back(muonTrack());
	getline(fin, temp, ',' );
	muonTracks[index].time = atof(temp.c_str());


	while (fin.good())
	{
		//muonTracks.push_back(muonTrack());
		//getline(fin, dummy, ',' );
		//muonTracks[index].time = atof(dummy.c_str());
        index++;
		getline(fin, temp, ',' );
		muonTracks[index].inDirX = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].inDirY = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].inDirZ = atof(temp.c_str());

		getline(fin, temp, ',' );
		muonTracks[index].inPointX = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].inPointY = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].inPointZ = atof(temp.c_str());

		getline(fin, temp, ',' );
		muonTracks[index].outDirX = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].outDirY = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].outDirZ = atof(temp.c_str());

		getline(fin, temp, ',' );
		muonTracks[index].outPointX = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].outPointY = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].outPointZ = atof(temp.c_str());

		getline(fin, temp, ',' );
		muonTracks[index].momentumIn = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].momentumOut = atof(temp.c_str());

		getline(fin, temp, ',' );
		muonTracks[index].scatAngleX = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].scatAngleY = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].scatAngleTotal = atof(temp.c_str());

		// POCA x y and z
		getline(fin, temp, ',' );
		muonTracks[index].pocaX = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].pocaY = atof(temp.c_str());
		getline(fin, temp, ',' );
		muonTracks[index].pocaZ = atof(temp.c_str());

		muonTracks.push_back(muonTrack());
		getline(fin, temp, ',' );
		muonTracks[index].time = atof(temp.c_str());

	}
	return index;
}

void initializeLookups(int * lowScatlookupX, int * lowScatlookupY, int * lowScatlookupZ, int length)
{
	for(int i = 0; i < length; i++)
	{
		lowScatlookupX[i] = 0;
		lowScatlookupY[i] = 0;
		lowScatlookupZ[i] = 0;
	}
}
