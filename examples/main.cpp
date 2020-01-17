/*

                          3D Zernike Moments
    Copyright (C) 2003 by Computer Graphics Group, University of Bonn
           http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/

Code by Marcin Novotni:     marcin@cs.uni-bonn.de

for more information, see the paper:

@inproceedings{novotni-2003-3d,
    author = {M. Novotni and R. Klein},
    title = {3{D} {Z}ernike Descriptors for Content Based Shape Retrieval},
    booktitle = {The 8th ACM Symposium on Solid Modeling and Applications},
    pages = {216--225},
    year = {2003},
    month = {June},
    institution = {Universit\"{a}t Bonn},
    conference = {The 8th ACM Symposium on Solid Modeling and Applications, June 16-20, Seattle, WA}
}
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

/*
    This is a demonstration of how the classes may be used to generate the
    3D Zernike descriptors from a given input binary file containing the
    voxel grid representation of the object.

    Notice that in the below case, the file contains a cubic grid, i.e. the
    x-, y-, and z-dimensions are equal (such setting should satisfy most needs).
*/

#include <iostream>
#include <boost\filesystem.hpp>

// ---- local includes ----
#include "ZernikeDescriptor.h"
#include "BinvoxReader.h"


// reads a voxel grid from a binary file
template<class TIn, class TOut>
TOut* ReadGrid (const char* _fname, int& _dim)

{
    std::ifstream infile (_fname, std::ios_base::binary | std::ios_base::in);

    vector<double> tempGrid;
    TIn temp;

    // read the grid values
    while (infile.read ((char*)(&temp), sizeof (TIn)))
    {
        tempGrid.push_back ((double)temp);
    }

    int d = tempGrid.size ();
    double f = std::pow<double> ((double)d, 1.0/3.0);

    _dim = floor (f+0.5);

    TOut* result = new TOut [d];
    for (int i=0; i<d; ++i)
    {
        result[i] = tempGrid[i];
    }

    return result;
}



int main(int argc, char** argv)
{
	if (argc != 3)
	{
		std::cout << "Usage: ZernikeMoments " <<
			"directname " <<
			"<MaxOrder> " << std::endl;
		return 0;
	}
	boost::filesystem::path DirectName(argv[1]);
	std::vector<boost::filesystem::path> v;

	for (auto&& x : boost::filesystem::recursive_directory_iterator(DirectName)) {
		std::cout << x.path() << std::endl;
		auto pos = x.path().string().find_first_of(".");
		if ((x.path().string()[pos + 1] == 's') && (x.path().string()[pos + 2] == 'u') && (x.path().extension() == ".binvox")) {
			std::cout << x.path() << std::endl;
			v.push_back(x.path());
		}

	}

		//int d;
		//double* voxels = ReadGrid<float, double> (argv[1], d);

		// .inv file name
	for (size_t i=0; i<v.size()-1; ++i){

		std::string path{ v[i].string() }, invFName;

		vector<unsigned char> voxels;

		size_t dim{};

		if (!io::binvox::read_binvox(path, voxels, dim))
		{
			std::cerr << "Cannot read binvox" << std::endl;
			continue;
		}

		// Воксель содержащий значения из центра
		vector<double> double_voxels(voxels.size());

		for (size_t x{ 0 }; x < dim; x++)
		{
			for (size_t z{ 0 }; z < dim; z++)
			{
				for (size_t y{ 0 }; y < dim; y++)
				{
					double_voxels.at((x * dim + z) * dim + y) = voxels.at((x * dim + z) * dim + y);
				}
			}
		}

		auto pos = path.find_last_of(".");

		if (pos != std::string::npos && pos != 0)
		{
			invFName = path.substr(0, pos);
		}
		else
		{
			std::cerr << "No extension in input filename? : " << v[i].string() << std::endl;
		}

		// compute the zernike descriptors
		ZernikeDescriptor<double, double> zd(double_voxels.data(), dim, std::stoi(argv[2]));

		invFName += ".inv";

		std::string invFile = invFName;
		
		std::cout << "Saving invariants file: " << invFile << std::endl;

		// save them into an .inv file
		zd.SaveInvariants(invFile.c_str());
	}
}
