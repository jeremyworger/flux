/*
 * Flux: Converge towards a target image using a Monte Carlo algorithm.
 *		 Run "flux -h" for help.
 *
 * Copyright (C) 2013  Lester Hedges <lester.hedges+flux@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "lodepng.h"
#include "MersenneTwister.h"

using namespace std;

// FUNCTION PROTOTYPES
void printHelpMessage();
void parseCommandLineArguments(int, char**, ostringstream&, ostringstream&, long&, unsigned int&, double&, bool&, bool&, bool&);
void generateSamplePoints(vector <double>&, long&, unsigned int&, bool&);
void decodeImage(ostringstream&, vector <unsigned char>&, unsigned int&, unsigned int&);
void encodeImage(unsigned int&, ostringstream&, vector <unsigned char>&, unsigned int&, unsigned int&);
void swapPixels(vector <unsigned char>&, vector <vector <unsigned int> >&, unsigned int&, unsigned int&, MTRand&, unsigned int&, unsigned int&);
double computePixelValue(vector <unsigned char>&, vector <unsigned char>&, unsigned int&, unsigned int&);
void randomizeImage(vector <unsigned char>&, MTRand&, unsigned int&, unsigned int&);
void computeNeighbors(vector <vector <unsigned int> >&, unsigned int&, unsigned int&);

// BEGIN MAIN FUNCTION
int main (int argc, char** argv)
{
	// counters
	long steps = 0;						// trial iteration counter
	long nAccepted = 0;					// number of accepted trials
	unsigned int samples = 0;			// number of samples

	// file name buffers
	ostringstream targetFileName;		// file name of target image
	ostringstream startingFileName;		// file name of starting image

	// defaults
	long mcSweeps = 10000;				// total number of MC sweeps
	unsigned int frames = 1000;			// number of images to sample
	double temperature = 0.1;			// temperature of thermal bath
	bool isLogarithmic = false;			// whether sampling is performed at logarithmic intervals
	bool isReverse = false;				// whether running in reverse mode (dissolve image)
	bool isMonitor = false;				// whether to monitor acceptance statistics

	// absoulte pixel value of image block before and after trial
	double currentPixelValue, trialPixelValue, pixelValueChange;

	// width and height of image in pixels
	unsigned int width, height;

	// pixel indices
	unsigned int p1, p2;

	// Monte Carlo sweeps
	double sweeps = 0;

	// random number generator
	MTRand rng;

	// output directory name (none selected if empty)
	ostringstream directoryName;
	directoryName.str("");

	// parse command-line arguments
	parseCommandLineArguments(argc, argv, targetFileName, directoryName,
		mcSweeps, frames, temperature, isLogarithmic, isReverse, isMonitor);

	// file name variable for stat logging
	ostringstream fileName;

	// file pointer
	FILE *pFile;

	if (isMonitor)
	{
		fileName.str("");

		if (strcmp(directoryName.str().c_str(), "") == 0)
		{
			fileName << "log.txt";
		}
		else
		{
			fileName << directoryName.str() << "/log.txt";
		}

		// wipe existing log file
		pFile = fopen(fileName.str().c_str(), "w");
		fclose(pFile);
	}

	// array of sampling check points
	vector <double> samplePoints;

	// generate sampling points
	generateSamplePoints(samplePoints, mcSweeps, frames, isLogarithmic);

	// set inverse temperature
	double beta = 1.0/temperature;

	// initialize and read target image from file
	vector <unsigned char> targetImage, trialImage;

	// decode target image from file
	decodeImage(targetFileName, targetImage, width, height);

	// set total number of pixels in canvas
	unsigned int pixels = width*height;

	// generate neigbor list
	vector <vector <unsigned int> > neighbors(width*height, vector <unsigned int> (4));
	computeNeighbors(neighbors, width, height);

	// acceptance flag
	bool isAccepted;

	trialImage = targetImage;
	if (!isReverse) randomizeImage(trialImage, rng, width, height);

	// compute current distance from target
	double error = 0;

	for (unsigned int i=0;i<width*height;i++)
	{
		error += computePixelValue(trialImage, targetImage, i, i);
	}

	cout << "Starting image generation..." << endl;

	double tStep = 1.0/(width*height);

	// MAIN LOOP
	while (sweeps < mcSweeps)
	{
		steps++;
		sweeps += tStep;

		// reset acceptance flag
		isAccepted = false;

		// choose random pixel and neighbor
		p1 = rng.randInt(width*height - 1);
		p2 = neighbors[p1][rng.randInt(3)];

		// current difference between RGBA magnitude in trial and target images
		currentPixelValue = computePixelValue(trialImage, targetImage, p1, p2);

		// swap pixel RGBA values
		swapPixels(trialImage, neighbors, width, height, rng, p1, p2);

		// new difference between RGBA magnitude in trial and target images
		trialPixelValue = computePixelValue(trialImage, targetImage, p1, p2);

		// evaluate change in pixel value
		pixelValueChange = trialPixelValue - currentPixelValue;

		// reverse
		if (isReverse) pixelValueChange *= -1;

		// check whether move is accepted
		if (pixelValueChange < 0) isAccepted = true;
		else
		{
			if (rng() < exp(-beta*pixelValueChange)) isAccepted = true;
		}

		// copy pixels from trial to current image if accepted, do the opposite if rejected
		if (isAccepted)
		{
			nAccepted++;
			error += pixelValueChange;
		}
		else swapPixels(trialImage, neighbors, width, height, rng, p2, p1);

		// encode current image and write stats to stdout and file
		if (sweeps >= samplePoints[samples])
		{
			samples++;

			// write to stdout
			printf("Frame: %4d, error: %5.4f, acceptance: %5.4f\n", samples, error/pixels, ((double) nAccepted/steps));

			// write acceptance statistics to file
			if (isMonitor)
			{
				// write to log file
				pFile = fopen(fileName.str().c_str(), "a");
				fprintf(pFile, "%d %5.4f %5.4f\n", samples, error/pixels, ((double) nAccepted/steps));
				fclose(pFile);
			}

			// write current evolution
			encodeImage(samples, directoryName, trialImage, width, height);
		}
	}
	// END MAIN LOOP

	cout << "Complete!" << endl;

	return 0;
}
// END MAIN FUNCTION

// FUNCTION DEFINITIONS

// Print help message to stdout
void printHelpMessage()
{
	puts("Flux version 0.1.0, by Lester O. Hedges\n\n"
	     "Syntax: flux image [options]\n\n"

	     "Available options:\n"
	     " -h/--help                 : Print this help information\n"
	     " -i/--iterations <int>     : Number of iterations (MC sweeps)\n"
	     " -f/--frames <int>         : Number of frames to encode\n"
	     " -t/--temperature <double> : Temperature of thermal bath\n"
	     " -l/--log                  : Sample at logarithmic intervals\n"
	     " -r/--reverse              : Run in reverse mode (dissolve image)\n"
	     " -d/--directory            : Name of output directory (will be created)\n"
	     " -m/--monitor              : Monitor acceptance statistics\n");
}

// Parse arguments from command-line
void parseCommandLineArguments(int argc, char **argv, ostringstream &targetFileName, ostringstream &directoryName,
	long &sweeps, unsigned int &frames, double &temperature, bool &isLogarithmic, bool &isReverse, bool &isMonitor)
{
	int i = 1;
	string s;

	// check whether only arg is requesting help
	if (argc == 2)
	{
		if(strcmp(argv[1],"-h") == 0 || strcmp(argv[1],"--help") == 0)
		{
			printHelpMessage();
			exit(0);
		}
	}

	// must have at least 3 arguments
	if (argc < 2)
	{
		printHelpMessage();
		exit(1);
	}
	else
	{
		s = argv[i];
		targetFileName.str(s);
		i++;

		while (i < argc)
		{
			if(strcmp(argv[i],"-i") == 0 || strcmp(argv[i],"--iterations") == 0)
			{
				i++;
				sweeps = atoi(argv[i]);
			}

			else if(strcmp(argv[i],"-f") == 0 || strcmp(argv[i],"--frames") == 0)
			{
				i++;
				frames = atoi(argv[i]);
			}

			else if(strcmp(argv[i],"-d") == 0 || strcmp(argv[i],"--directory") == 0)
			{
				i++;
				s = argv[i];
				directoryName.str(s);

				struct stat st;
				string systemCommand;

				// check to see if folder already exists
				if (stat(directoryName.str().c_str(), &st) != 0)
				{
					// create top level directory for output
					systemCommand = "mkdir " + directoryName.str();
					system(systemCommand.c_str());

					cout << "Created output directory: " << directoryName.str() << "/\n" << endl;
				}
			}

			else if(strcmp(argv[i],"-t") == 0 || strcmp(argv[i],"--temperature") == 0)
			{
				i++;
				temperature = atof(argv[i]);
			}

			else if(strcmp(argv[i],"-l") == 0 || strcmp(argv[i],"--logarithmic") == 0)
			{
				isLogarithmic = true;
			}

			else if(strcmp(argv[i],"-r") == 0 || strcmp(argv[i],"--reverse") == 0)
			{
				isReverse = true;
			}

			else if(strcmp(argv[i],"-m") == 0 || strcmp(argv[i],"--monitor") == 0)
			{
				isMonitor = true;
			}

			else
			{
				cerr << "Error: unknown command-line parameter " << argv[i] << endl;
				cerr << "Run: for help run \"mc_charcoal -h\"" << endl;
				exit(1);
			}

			i++;
		}
	}
}

// Generate array of linear or logarithmic sampling points
void generateSamplePoints(vector <double> &samplePoints, long &sweeps, unsigned int &frames, bool &isLogarithmic)
{
	samplePoints.resize(frames);

	unsigned int i;
	if (isLogarithmic)
	{
		for (i=0;i<frames;i++)
		{
			samplePoints[i] = pow(sweeps, ((double) i / (frames-1)));
		}
	}
	else
	{
		double sampleInterval = sweeps/frames;

		for (i=0;i<frames;i++) samplePoints[i] = (i+1)*sampleInterval;
	}
}

// Decode png image from disk
void decodeImage(ostringstream &fileName, vector <unsigned char> &image, unsigned int &width, unsigned int &height)
{
	// decode
	unsigned error = lodepng::decode(image, width, height, fileName.str().c_str());

	//if there's an error, display it
	if(error) cout << "Decoder error " << error << ": " << lodepng_error_text(error) << endl;

	// N.B. lodepng reads in 4 bytes per pixel, i.e. data is stored as RGBARGBA...
}

// Encode png image to disk
void encodeImage(unsigned int &frame, ostringstream &directoryName, vector <unsigned char> &image, unsigned int &width, unsigned int &height)
{
	ostringstream fileName,num;
	num.str("");
	num.width(5);
	num.fill('0');
	num << right << frame;

	fileName.str("");

	if (strcmp(directoryName.str().c_str(), "") == 0)
	{
		fileName << "Frame_" << num.str() << ".png";
	}
	else fileName << directoryName.str() << "/Frame_" << num.str() << ".png";

	// encode the image
	unsigned error = lodepng::encode(fileName.str().c_str(), image, width, height);

	// if there's an error, display it
	if(error) cout << "Encoder error " << error << ": "<< lodepng_error_text(error) << endl;
}

// Swap two neighboring pixels
void swapPixels(vector <unsigned char> &trialImage, vector <vector <unsigned int> > &neighbors,
	unsigned int &width, unsigned int &height, MTRand &rng, unsigned int &p1, unsigned int &p2)
{

	unsigned char tmp;

	// swap pixel RGBA values
	for (unsigned int i=0;i<4;i++)
	{
		tmp = trialImage[4 * p1 + i];
		trialImage[4 * p1 + i] = trialImage[4 * p2 + i];
		trialImage[4 * p2 + i] = tmp;
	}
}

// Compute difference in RGBA pixel magnitude between two images
double computePixelValue(vector <unsigned char> &image, vector <unsigned char> &target, unsigned int &p1, unsigned int &p2)
{
	double m1 = 0;
	double m2 = 0;
	double norm = 255*255;

	for (unsigned int i=0;i<4;i++)
	{
		m1 += (int(image[4 * p1 + i]) - int(target[4 * p1 + i])) * (int(image[4 * p1 + i]) - int(target[4 * p1 + i]))/norm;
		m2 += (int(image[4 * p2 + i]) - int(target[4 * p2 + i])) * (int(image[4 * p2 + i]) - int(target[4 * p2 + i]))/norm;
	}

	return (sqrt(m1) + sqrt(m2))/4.0;
}

// Randomize pixels
void randomizeImage(vector <unsigned char> &image, MTRand &rng, unsigned int &width, unsigned int &height)
{
	unsigned int i,j;
	unsigned int pixel;
	unsigned int nPixels = width*height;
	vector <unsigned int> pixels(nPixels);

	// copy original image
	vector <unsigned char> originalImage = image;

	// fill pixel index array
	for (i=0;i<nPixels;i++) pixels[i] = i;

	for (i=0;i<nPixels;i++)
	{
		pixel = pixels[rng.randInt(nPixels-i-1)];
		pixels[pixel] = pixels[nPixels-i-1];

		// copy pixel into position
		for (j=0;j<4;j++)
		{
			image[4 * i + j] = originalImage[4 * pixel + j];
		}
	}
}

// Compute nearest neighbors of each pixel using periodic boundaries
void computeNeighbors(vector <vector <unsigned int> > &neighbors, unsigned int &width, unsigned int &height)
{
	int x, y;
	unsigned int nPixels = width*height;

	div_t coord;

	for (unsigned int i=0;i<nPixels;i++)
	{
		coord = div(int(i), int(width));

		x = coord.rem;
		y = coord.quot;

		neighbors[i][0] = width*y + (x-1+width)%width;
		neighbors[i][1] = width*y + (x+1)%width;
		neighbors[i][2] = ((y-1+height)%height)*width + x;
		neighbors[i][3] = ((y+1)%height)*width + x;
	}
}
