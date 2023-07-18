#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <stdbool.h>

typedef struct trajectory
{
	int atomID, atomType, molType, ix, iy, iz;
	float x, y, z;
	int isEndGroup;
} TRAJECTORY;

typedef struct vector
{
	float x1, y1, z1;
	float x2, y2, z2;
	float xc, yc, zc;
} VECTOR;

typedef struct datafileInfo
{
	int nAtoms, nBonds, nAngles, nDihedrals, nImpropers;
	int nAtomTypes, nBondTypes, nAngleTypes, nDihedralTypes, nImproperTypes;
} DATAFILE_INFO;

typedef struct datafile_atoms
{
	int resNumber;
	char resName[6], atomName[6], atomType2[6], molName[6];

	int id, molType, atomType;
	float charge, x, y, z;
} DATA_ATOMS;

typedef struct datafile_bonds
{
	int id, bondType, atom1, atom2;
} DATA_BONDS;

typedef struct datafile_angles
{
	int id, angleType, atom1, atom2, atom3;
} DATA_ANGLES;

typedef struct datafile_dihedrals
{
	int id, dihedralType, atom1, atom2, atom3, atom4;
} DATA_DIHEDRALS;

typedef struct datafile_impropers
{
	int id, improperType, atom1, atom2, atom3, atom4;
} DATA_IMPROPERS;

typedef struct simulationBoundary
{
	float xlo, xhi, ylo, yhi, zlo, zhi;
	float xLength, yLength, zLength;
} SIMULATION_BOUNDARY;

typedef struct rdf
{
	float rlo, rhi, gofr;
} RDF;

typedef struct stats
{
	float average, standardDeviation;
} STATS;

typedef struct orderParameterBins
{
	float orderParameter, rlo, rhi, count;
} ORDERPARAMETER_BINS;

typedef struct distanceBins
{
	float rlo, rhi, count;
} DIST_BINS;

typedef struct angleDistribution
{
	float thetaLo, thetaHi, count;
} ANGLE_DISTRIBUTION;

SIMULATION_BOUNDARY readDumpBoundary (FILE *file_dump, SIMULATION_BOUNDARY boundary)
{
	rewind (file_dump);
	char lineString[2000];

	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 2000, file_dump);
	}

	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.xlo, &boundary.xhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.ylo, &boundary.yhi);
	fgets (lineString, 2000, file_dump);
	sscanf (lineString, "%f %f\n", &boundary.zlo, &boundary.zhi);
	rewind (file_dump);

	boundary.xLength = boundary.xhi - boundary.xlo;
	boundary.yLength = boundary.yhi - boundary.ylo;
	boundary.zLength = boundary.zhi - boundary.zlo;

	printf("xlo xhi %f %f\nylo yhi %f %f\nzlo zhi %f %f\n\n", boundary.xlo, boundary.xhi, boundary.ylo, boundary.yhi, boundary.zlo, boundary.zhi);

	return boundary;
}

TRAJECTORY *readTimestep (FILE *file_dump, TRAJECTORY *atoms, int nAtomEntries, SIMULATION_BOUNDARY *boundary)
{
	char lineString[2000];
	int currentAtomID = 1;

	for (int i = 0; i < 5; ++i) {
		fgets (lineString, 2000, file_dump); }

	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).xlo, &(*boundary).xhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).ylo, &(*boundary).yhi);
	fgets (lineString, 2000, file_dump); sscanf (lineString, "%f %f\n", &(*boundary).zlo, &(*boundary).zhi);
	fgets (lineString, 2000, file_dump);

	for (int i = 0; i < nAtomEntries; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);
		sscanf (lineString, "%d %d %f %f %f %d %d %d\n", &atoms[currentAtomID - 1].atomID, &atoms[currentAtomID - 1].atomType, &atoms[currentAtomID - 1].x, &atoms[currentAtomID - 1].y, &atoms[currentAtomID - 1].z, &atoms[currentAtomID - 1].ix, &atoms[currentAtomID - 1].iy, &atoms[currentAtomID - 1].iz);
		atoms[currentAtomID - 1].isEndGroup = 0;
	}

	return atoms;
}

int **initNeighIDs (int **neighborIDs, int nAtoms, int coordinationNumber)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		for (int j = 0; j < coordinationNumber; ++j)
		{
			neighborIDs[i][j] = -1;
		}
	}

	return neighborIDs;
}

TRAJECTORY *initializeAtoms (TRAJECTORY *atoms, int nAtoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		atoms[i].atomID = 0;
		atoms[i].atomType = 0;
		atoms[i].molType = 0;
		atoms[i].ix = 0;
		atoms[i].iy = 0;
		atoms[i].iz = 0;
		atoms[i].x = 0;
		atoms[i].y = 0;
		atoms[i].z = 0;
		atoms[i].isEndGroup = 0;
	}

	return atoms;
}

float translatePeriodicDistance (float x1, float x2, float simBoxLength, float newR)
{
	if (fabs (x1 - x2) > (simBoxLength / (float)2))
	{
		if (x1 >= x2) {
			newR = x1 - simBoxLength; }
		else if (x1 < x2) {
			newR = x1 + simBoxLength; }

		return newR;
	}
	else
	{
		return x1;
	}
}

int countNAtoms (int *nAtomEntries, const char *inputFileName)
{
	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));
	FILE *file_dump;

	if (strstr (inputFileName, ".xz")) {
		snprintf (pipeString, 500, "xzcat %s", inputFileName);
		file_dump = popen (pipeString, "r"); }
	else {
		file_dump = fopen (inputFileName, "r"); }

	int nAtoms, currentAtomID, nAtomsFixed;
	char lineString[2000];

	for (int i = 0; i < 4; ++i) {
		fgets (lineString, 2000, file_dump); }

	sscanf (lineString, "%d\n", &nAtoms);
	(*nAtomEntries) = nAtoms;
	
	if (strstr (inputFileName, ".xz")) {
		pclose (file_dump);
		file_dump = popen (pipeString, "r"); }
	else {
		rewind (file_dump); }

	nAtomsFixed = nAtoms;

	for (int i = 0; i < 9; ++i) {
		fgets (lineString, 2000, file_dump); }

	for (int i = 0; i < nAtoms; ++i)
	{
		fgets (lineString, 2000, file_dump);
		sscanf (lineString, "%d\n", &currentAtomID);

		if (currentAtomID > nAtoms) {
			nAtomsFixed = currentAtomID; }
	}

	printf("Number of atom entries in the dump file: %d\nTotal number of atoms present in the simulation: %d\n\n", nAtoms, nAtomsFixed);

	if (strstr (inputFileName, ".xz")) {
		pclose (file_dump); }
	else {
		fclose (file_dump); }

	return nAtomsFixed;
}

typedef struct pairPairDistance
{
	float distance;
	int atom1, atom2, index;
} PAIR_PAIR_DISTANCE;

int **findNeighbors2 (int **neighborIDs, PAIR_PAIR_DISTANCE *pairInfo, int i, int coordinationNumber, int nAtoms, int atomType1, int atomType2, TRAJECTORY *atoms)
{
	PAIR_PAIR_DISTANCE *nearestAtomDistances, max;
	nearestAtomDistances = (PAIR_PAIR_DISTANCE *) malloc (coordinationNumber * sizeof (PAIR_PAIR_DISTANCE));

	for (int j = 0; j < coordinationNumber; ++j)
	{
		nearestAtomDistances[j].distance = pairInfo[j].distance;
		nearestAtomDistances[j].atom1 = pairInfo[j].atom1;
		nearestAtomDistances[j].atom2 = pairInfo[j].atom2;
	}

	for (int l = 0; l < nAtoms; ++l)
	{
		if ((atoms[l].atomType == atomType1) || (atomType1 == -1))
		{
				for (int j = 0; j < coordinationNumber; ++j)
				{
					if (j == 0) {
						max.distance = nearestAtomDistances[j].distance;
						max.atom1 = nearestAtomDistances[j].atom1;
						max.atom2 = nearestAtomDistances[j].atom2;
						max.index = j; }
					else
					{
						if (nearestAtomDistances[j].distance > max.distance) {
							max.distance = nearestAtomDistances[j].distance; 
							max.atom1 = nearestAtomDistances[j].atom1;
							max.atom2 = nearestAtomDistances[j].atom2;
							max.index = j; }
					}
				}

				if (pairInfo[l].distance < max.distance) {
					nearestAtomDistances[max.index].distance = pairInfo[l].distance;
					nearestAtomDistances[max.index].atom1 = pairInfo[l].atom1;
					nearestAtomDistances[max.index].atom2 = pairInfo[l].atom2; }
		}
	}

/*	printf("\n");
*/	for (int m = 0; m < coordinationNumber; ++m) {
		neighborIDs[i][m] = (int)nearestAtomDistances[m].atom2;
/*		printf("%d %d\n", i + 1, nearestAtomDistances[m].atom2);
		usleep (100000); 
*/	}

	return neighborIDs;
}

float computeDistance (TRAJECTORY *atoms, int i, int j, SIMULATION_BOUNDARY boundary)
{
	float xLength = (boundary.xhi - boundary.xlo), yLength = (boundary.yhi - boundary.ylo), zLength = (boundary.zhi - boundary.zlo);
	float distance;
	float newX, newY, newZ;

	newX = translatePeriodicDistance (atoms[i].x, atoms[j].x, xLength, newX);
	newY = translatePeriodicDistance (atoms[i].y, atoms[j].y, yLength, newY);
	newZ = translatePeriodicDistance (atoms[i].z, atoms[j].z, zLength, newZ);

	distance = sqrt (
		(newX - atoms[j].x) * (newX - atoms[j].x) +
		(newY - atoms[j].y) * (newY - atoms[j].y) +
		(newZ - atoms[j].z) * (newZ - atoms[j].z)
		);

	return distance;
}

int **findNeighbors (TRAJECTORY *atoms, int **neighborIDs, int nAtoms, SIMULATION_BOUNDARY boundary, int coordinationNumber, int atomType1, int atomType2)
{
	PAIR_PAIR_DISTANCE *pairInfo;
	pairInfo = (PAIR_PAIR_DISTANCE *) malloc (nAtoms * sizeof (PAIR_PAIR_DISTANCE));

	for (int i = 0; i < nAtoms; ++i)
	{
		if ((atoms[i].atomType == atomType1) || (atomType1 == -1))
		{
			for (int j = 0; j < nAtoms; ++j)
			{
				if (i != j)
				{
					if ((atoms[j].atomType == atomType2) || (atomType2 == -1))
					{
						pairInfo[j].distance = computeDistance (atoms, i, j, boundary);
						pairInfo[j].atom1 = i + 1;
						pairInfo[j].atom2 = j + 1;
					}
					else
					{
						pairInfo[j].distance = 9999;
						pairInfo[j].atom1 = i + 1;
						pairInfo[j].atom2 = j + 1;
					}
				}
				if (i == j)
				{
					pairInfo[j].distance = 9999;
					pairInfo[j].atom1 = i + 1;
					pairInfo[j].atom2 = j + 1;
				}
			}

			neighborIDs = findNeighbors2 (neighborIDs, pairInfo, i, coordinationNumber, nAtoms, atomType1, atomType2, atoms);
		}

	}

	return neighborIDs;
}

float *computeNonBondedOrderParameter (float *nonBondedAngle, TRAJECTORY *atoms, int nAtoms, int **neighborIDs, int atomType1, int atomType2, int coordinationNumber)
{
	float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4;
	float dotProduct, magnitude1, magnitude2, cosTheta, theta;
	int denominator;

	for (int i = 0; i < nAtoms; ++i)
	{
		theta = 0;
		denominator = 0;

		if ((atoms[i].atomType == atomType1) || (atomType1 == -1))
		{
			for (int j = 0; j < coordinationNumber; ++j)
			{
				for (int k = (j + 1); k < coordinationNumber; ++k)
				{
					x1 = atoms[i].x;
					y1 = atoms[i].y;
					z1 = atoms[i].z;

					x2 = atoms[neighborIDs[i][j]].x;
					y2 = atoms[neighborIDs[i][j]].y;
					z2 = atoms[neighborIDs[i][j]].z;

					x3 = atoms[i].x;
					y3 = atoms[i].y;
					z3 = atoms[i].z;

					x4 = atoms[neighborIDs[i][k]].x;
					y4 = atoms[neighborIDs[i][k]].y;
					z4 = atoms[neighborIDs[i][k]].z;

					dotProduct = ((x2 - x1) * (x4 - x3)) + ((y2 - y1) * (y4 - y3)) + ((z2 - z1) * (z4 - z3)); 
					magnitude1 = ((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)) + ((z2 - z1) * (z2 - z1)); 
					magnitude2 = ((x4 - x3) * (x4 - x3)) + ((y4 - y3) * (y4 - y3)) + ((z4 - z3) * (z4 - z3)); 

					cosTheta = dotProduct / (sqrt (magnitude1) * sqrt (magnitude2)); 
					theta += acosf (cosTheta); 
					denominator++;
				}
			}

			nonBondedAngle[i] = theta / (float)denominator;
		}
		else
		{
			nonBondedAngle[i] = -9999;
		}

	}

	return nonBondedAngle;
}

float *initNonBondedAngleAverage (float *nonBondedAngleAverage, int nAtoms, int atomType1, TRAJECTORY *atoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
			nonBondedAngleAverage[i] = 0;
	}

	return nonBondedAngleAverage;
}

float *sumNonBondedAngles (float *nonBondedAngleAverage, float *nonBondedAngle, int nAtoms, int atomType1, int atomType2, TRAJECTORY *atoms)
{
	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == atomType1)
		{
			nonBondedAngleAverage[i] += nonBondedAngle[i];
		}
	}

	return nonBondedAngleAverage;
}

void printDump (float *nonBondedAngle, int nAtoms, TRAJECTORY *atoms, int atomType1, int atomType2, FILE *file_nonbonded_dump, int currentTimestep)
{
	fprintf(file_nonbonded_dump, "Timestep: %d\nNumber of atoms: %d\nAtom type 1: %d\nAtom type 2: %d\nEntries: AtomID angle(degrees)\n", currentTimestep, nAtoms, atomType1, atomType2);

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == atomType1)
		{
			fprintf(file_nonbonded_dump, "%d %f\n", (i + 1), (nonBondedAngle[i] * 57.29578));
		}
	}
}

void printFinalAverage (float *nonBondedAngleAverage, int nAtoms, TRAJECTORY *atoms, int atomType1, int atomType2, FILE *file_nonbonded_final, int currentTimestep)
{
	fprintf(file_nonbonded_final, "#This file contains average angles\nTimestep: %d\nNumber of atoms: %d\nAtom type 1: %d\nAtom type 2: %d\nEntries: AtomID angle(degrees)\n", currentTimestep, nAtoms, atomType1, atomType2);

	for (int i = 0; i < nAtoms; ++i)
	{
		if (atoms[i].atomType == atomType1)
		{
			fprintf(file_nonbonded_final, "%d %f\n", (i + 1), (nonBondedAngleAverage[i] * 57.29578));
		}
	}
}

ANGLE_DISTRIBUTION *initAngleDistribution (ANGLE_DISTRIBUTION *angleDistribution, int angle_nBins, float angle_binWidth)
{
	for (int i = 0; i < angle_nBins; ++i)
	{
		if (i == 0)
		{
			angleDistribution[i].thetaLo = 0;
			angleDistribution[i].thetaHi = angle_binWidth;
			angleDistribution[i].count = 0;
		}
		else
		{
			angleDistribution[i].thetaLo = angleDistribution[i - 1].thetaHi;
			angleDistribution[i].thetaHi = angleDistribution[i].thetaLo + angle_binWidth;
			angleDistribution[i].count = 0;
		}
	}

	return angleDistribution;
}

ANGLE_DISTRIBUTION *computeNonBondedAngleDistribution (ANGLE_DISTRIBUTION *nonBondedAngleDistribution, float *nonBondedAngle, int nAtoms, float angle_binWidth, int angle_nBins)
{
	float thetaDegrees;

	for (int i = 0; i < nAtoms; ++i)
	{
		if (nonBondedAngle[i] > -9000)
		{
			for (int j = 0; j < angle_nBins; ++j)
			{
				thetaDegrees = nonBondedAngle[i] * 57.29578;

				if ((thetaDegrees <= nonBondedAngleDistribution[j].thetaHi) && (thetaDegrees > nonBondedAngleDistribution[j].thetaLo))
				{
					nonBondedAngleDistribution[j].count += 1;
				}
			}
		}
	}
	return nonBondedAngleDistribution;
}

ANGLE_DISTRIBUTION *computeNonBondedAngleDistributionAverage (ANGLE_DISTRIBUTION *nonBondedAngleDistributionAverage, ANGLE_DISTRIBUTION *nonBondedAngleDistribution, float angle_binWidth, int angle_nBins)
{
	for (int i = 0; i < angle_nBins; ++i) {
		nonBondedAngleDistributionAverage[i].count += nonBondedAngleDistribution[i].count; }

	return nonBondedAngleDistributionAverage;
}

void printAngleDistribution (ANGLE_DISTRIBUTION *nonBondedAngleDistribution, int angle_nBins, FILE *file_dist_rt, int currentTimestep)
{
	fprintf(file_dist_rt, "Timestep: %d\nEntries: thetaLo thetaHi count\n", currentTimestep);

	for (int i = 0; i < angle_nBins; ++i)
	{
		fprintf(file_dist_rt, "%f %f %.0f\n", nonBondedAngleDistribution[i].thetaLo, nonBondedAngleDistribution[i].thetaHi, nonBondedAngleDistribution[i].count);
	}
}

void printAngleDistributionAverage (ANGLE_DISTRIBUTION *nonBondedAngleDistributionAverage, int angle_nBins, FILE *file_dist, int currentTimestep)
{
	fprintf(file_dist, "Timestep: %d\nEntries: thetaLo thetaHi count\n", currentTimestep);

	for (int i = 0; i < angle_nBins; ++i)
	{
		fprintf(file_dist, "%f %f %.0f\n", nonBondedAngleDistributionAverage[i].thetaLo, nonBondedAngleDistributionAverage[i].thetaHi, nonBondedAngleDistributionAverage[i].count);
	}
}

int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		fprintf(stdout, "REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = program\n {~} argv[1] = input dump file name\n {~} argv[2] = atom type 1\n {~} argv[3] = atom type 2\n {~} argv[4] = Number of neighbors to consider\n {~} argv[5] = distance bin width\n {~} argv[6] = angle bin width.\n\n");
		fflush (stdout);
		exit (1);
	}

	FILE *file_dump, *file_dist, *file_dist_rt, *file_nonbonded_dump, *file_nonbonded_final;
	char *pipeString;
	pipeString = (char *) malloc (500 * sizeof (char));

	if (strstr (argv[1], ".xz")) {
		snprintf (pipeString, 500, "xzcat %s", argv[1]);
		file_dump = popen (pipeString, "r"); }
	else {
		file_dump = fopen (argv[1], "r"); }

	file_dist = fopen ("orderParameter.nonBonded.distribution", "w");
	file_nonbonded_dump = fopen ("orderParameter.nonBonded.dump", "w");
	file_nonbonded_final = fopen ("orderParameter.nonBonded.final", "w");
	file_dist_rt = fopen ("orderParameter.nonBonded.distribution.rt", "w");

	SIMULATION_BOUNDARY boundary;
	boundary = readDumpBoundary (file_dump, boundary);

	int nAtomEntries, nAtoms = countNAtoms (&nAtomEntries, argv[1]), atomType1 = atoi (argv[2]), atomType2 = atoi (argv[3]), file_status;

	if (strstr (argv[1], ".xz")) {
		pclose (file_dump);
		file_dump = popen (pipeString, "r"); }
	else {
		rewind (file_dump); }

	float angle_binWidth = atof (argv[6]), dist_binWidth = atof (argv[5]);
	int dist_nBins = ceil (30.0 / dist_binWidth), angle_nBins = ceil (360.0 / angle_binWidth);
	int coordinationNumber = atoi (argv[4]);

	TRAJECTORY *atoms;
	atoms = (TRAJECTORY *) malloc (nAtoms * sizeof (TRAJECTORY));
	printf("Number of atoms in the trajectory file: %d\n", nAtoms);

	int **neighborIDs;
	neighborIDs = (int **) malloc (nAtoms * sizeof (int *));

	float *nonBondedAngle, *nonBondedAngleAverage;
	nonBondedAngle = (float *) malloc (nAtoms * sizeof (float));
	nonBondedAngleAverage = (float *) malloc (nAtoms * sizeof (float));
	nonBondedAngleAverage = initNonBondedAngleAverage (nonBondedAngleAverage, nAtoms, atomType1, atoms);

	ANGLE_DISTRIBUTION *nonBondedAngleDistribution, *nonBondedAngleDistributionAverage;
	nonBondedAngleDistribution = (ANGLE_DISTRIBUTION *) malloc (angle_nBins * sizeof (ANGLE_DISTRIBUTION));
	nonBondedAngleDistributionAverage = (ANGLE_DISTRIBUTION *) malloc (angle_nBins * sizeof (ANGLE_DISTRIBUTION));
	nonBondedAngleDistribution = initAngleDistribution (nonBondedAngleDistribution, angle_nBins, angle_binWidth);
	nonBondedAngleDistributionAverage = initAngleDistribution (nonBondedAngleDistributionAverage, angle_nBins, angle_binWidth);

	for (int i = 0; i < nAtoms; ++i) {
		neighborIDs[i] = (int *) malloc (coordinationNumber * sizeof (int)); }

	atoms = initializeAtoms (atoms, nAtoms);

	if (strstr (argv[1], ".xz")) {
		pclose (file_dump);
		file_dump = popen (pipeString, "r"); }
	else {
		rewind (file_dump); }

	file_status = fgetc (file_dump);

	int currentTimestep = 0;

	while (file_status != EOF)
	{
		fprintf(stdout, "Scanning timestep: %d                           \r", currentTimestep);
		fflush (stdout);

		neighborIDs = initNeighIDs (neighborIDs, nAtoms, coordinationNumber);
		nonBondedAngleDistribution = initAngleDistribution (nonBondedAngleDistribution, angle_nBins, angle_binWidth);

		atoms = readTimestep (file_dump, atoms, nAtomEntries, &boundary);
		neighborIDs = findNeighbors (atoms, neighborIDs, nAtoms, boundary, coordinationNumber, atomType1, atomType2);
		nonBondedAngle = computeNonBondedOrderParameter (nonBondedAngle, atoms, nAtoms, neighborIDs, atomType1, atomType2, coordinationNumber);
		nonBondedAngleDistribution = computeNonBondedAngleDistribution (nonBondedAngleDistribution, nonBondedAngle, nAtoms, angle_binWidth, angle_nBins);
		nonBondedAngleDistributionAverage = computeNonBondedAngleDistributionAverage (nonBondedAngleDistributionAverage, nonBondedAngleDistribution, angle_binWidth, angle_nBins);

		printDump (nonBondedAngle, nAtoms, atoms, atomType1, atomType2, file_nonbonded_dump, currentTimestep);
		printAngleDistribution (nonBondedAngleDistribution, angle_nBins, file_dist_rt, currentTimestep);

		nonBondedAngleAverage = sumNonBondedAngles (nonBondedAngleAverage, nonBondedAngle, nAtoms, atomType1, atomType2, atoms);

		currentTimestep++;
		file_status = fgetc (file_dump);
	}

	printFinalAverage (nonBondedAngleAverage, nAtoms, atoms, atomType1, atomType2, file_nonbonded_final, currentTimestep);
	printAngleDistributionAverage (nonBondedAngleDistributionAverage, angle_nBins, file_dist, currentTimestep);

	fclose (file_dump);
	fclose (file_dist);
	fclose (file_dist_rt);
	fclose (file_nonbonded_dump);
	fclose (file_nonbonded_final);
	return 0;
}