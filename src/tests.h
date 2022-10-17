/*
 * tests.h
 *
 *  Created on: Apr 26, 2021
 *      Author: frederik
 */

#ifndef SRC_TESTS_H_
#define SRC_TESTS_H_


#include "includes_and_names.h"
#include <random>
#include "Array3D.h"

void testModulo()
{
	cout << 1 % 4 << endl;
	cout << 2 % 4 << endl;
	cout << 3 % 4 << endl;
	cout << 4 % 4 << endl;
	cout << 5 % 4 << endl;
	cout << 16 % 4 << endl;
	cout << 17 % 4 << endl;
	// Doesn't work:
	//cout << 10%0;
}

void testModuloWithDouble()
{
	cout << fmod(10.5, 2.5) << endl;
}

void testOutputFileStream()
{
	ofstream outFile("testFolder/testfile.csv");
	if (outFile)
	{
		cout << "File is open." << endl;
		outFile << "abc,1,2,3" << endl << "q,w,e,r";
		outFile.close();
	}
	else
		cout << "File wasn't opened." << endl;
}

void fillArrayRandomNumbers(Array3D_d& array, double minVal, double maxVal)
{
	for (uint i{0}; i<array.getSize(); ++i)
	{
		array(i) = (double)rand() / (double)RAND_MAX * (maxVal - minVal) + minVal ;
	}
}

void testArrayTraversalOneOrThreeIndex()
{
	Array3D_d A1(200, 200, 200);
	Array3D_d A2(200, 200, 200);
	Array3D_d A3(200, 200, 200);
	srand(time(NULL));
	fillArrayRandomNumbers(A1, 0, 10);
	fillArrayRandomNumbers(A2, 0, 10);
	fillArrayRandomNumbers(A3, 0, 10);

	Clock timer;
	for (uint i{0}; i<A1.getLength(); ++i)
		for (uint j{0}; j<A1.getWidth(); ++j)
			for (uint k{0}; k<A1.getHeight(); ++k)
				A3(i,j,k) = A1(i,j,k) + A2(i,j,k);
	Time t1 = timer.restart();
	for (uint i{0}; i<A1.getSize(); ++i)
		A3(i) = A1(i) + A2(i);
	Time t2 = timer.getElapsedTime();

	cout << "Triple index: " << t1.asSeconds() << " s" << endl;
	cout << "Single index: " << t2.asSeconds() << " s" << endl;
}

void testVectorAllocationSpeed()
{
	Clock clock;
	vector<double> v1(10);
	v1[9] = 2.3;
	Time t1 = clock.restart();
	vector<double> v2(10000000);
	v2[9999999] = 7.7;
	Time t2 = clock.getElapsedTime();

	cout << "t1: " << t1.asMicroseconds() << endl;
	cout << "t2: " << t2.asMicroseconds() << endl;
}

void testVectorAllocationCapacity()
{
	vector<double> v1(100);
	vector<double> v2(123);
	vector<double> v3(337);

	cout << "v1: size = " << v1.size() << ", capacity = " << v1.capacity() << endl;
	cout << "v2: size = " << v2.size() << ", capacity = " << v2.capacity() << endl;
	cout << "v3: size = " << v3.size() << ", capacity = " << v3.capacity() << endl;
	cout << "Now adding one element to each vector..." << endl;
	v1.push_back(1.);
	v2.push_back(1.);
	v3.push_back(1.);
	cout << "v1: size = " << v1.size() << ", capacity = " << v1.capacity() << endl;
	cout << "v2: size = " << v2.size() << ", capacity = " << v2.capacity() << endl;
	cout << "v3: size = " << v3.size() << ", capacity = " << v3.capacity() << endl;
}

void testDivision()
{
	cout << "2/3 = " << 2/3 << endl;
	cout << "2.0/3 = " << 2.0/3 << endl;
	cout << "2./3 = " << 2./3 << endl;
	cout << "2%3 = " << 2%3 << endl;
}

void fillArrayRandomNumbers(vector<double>& array, double minVal, double maxVal)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(minVal, maxVal);
	for (double& d : array)
	{
		d = distribution(generator);
	}
}

void testPerformanceArrayLoop()
{
	uint N = 10000000;
	vector<double> rho(N);
	vector<double> rho_u(N);
	vector<double> rho_v(N);
	vector<double> rho_w(N);
	vector<double> E(N);

	vector<double> u(N);
	vector<double> v(N);
	vector<double> w(N);
	vector<double> p(N);
	vector<double> T(N);
	vector<double> mu(N);
	vector<double> kappa(N);

	fillArrayRandomNumbers(rho  , 0, 10);
	fillArrayRandomNumbers(rho_u, 0, 10);
	fillArrayRandomNumbers(rho_v, 0, 10);
	fillArrayRandomNumbers(rho_w, 0, 10);
	fillArrayRandomNumbers(E    , 0, 10);

	bool allInOneLoop = false;
	if (allInOneLoop)
	{
		Clock clock;
		for (uint i{0}; i<N; ++i)
		{
			u[i] = rho_u[i] / rho[i];
			v[i] = rho_v[i] / rho[i];
			w[i] = rho_w[i] / rho[i];
			p[i] = 0.4 * ( E[i] - rho[i]/2 * ( u[i]*u[i] + v[i]*v[i] + w[i]*w[i] ));
			T[i] = ( 1.4*p[i] - rho[i] ) / ( 1+rho[i] );
			mu[i] = pow( 1 + T[i], 1.5 ) *1.3 / ( 200*( T[i]+1.3 ) );
			kappa[i] = mu[i] / ( 0.4*0.72 );
		}
		Time time = clock.getElapsedTime();
		cout << "All in one loop. Time: " << time.asSeconds() << " s" << endl;
	}
	else
	{
		Clock clock;
		for (uint i{0}; i<N; ++i)
			u[i] = rho_u[i] / rho[i];
		for (uint i{0}; i<N; ++i)
			v[i] = rho_v[i] / rho[i];
		for (uint i{0}; i<N; ++i)
			w[i] = rho_w[i] / rho[i];
		for (uint i{0}; i<N; ++i)
			p[i] = 0.4 * ( E[i] - rho[i]/2 * ( u[i]*u[i] + v[i]*v[i] + w[i]*w[i] ));
		for (uint i{0}; i<N; ++i)
			T[i] = ( 1.4*p[i] - rho[i] ) / ( 1+rho[i] );
		for (uint i{0}; i<N; ++i)
			mu[i] = pow( 1 + T[i], 1.5 ) *1.3 / ( 200*( T[i]+1.3 ) );
		for (uint i{0}; i<N; ++i)
			kappa[i] = mu[i] / ( 0.4*0.72 );

		Time time = clock.getElapsedTime();
		cout << "Multiple loops. Time: " << time.asSeconds() << " s" << endl;
	}
}

void testArray3D()
{
	Array3D_d A(3,2,1);
	A(0,0,0) = 7.7;
	A(0,0,5) = A(0,0,0) + 1.3;
	cout << A(0,0,5);
}

void testArrayWithSettings()
{
	ConfigSettings settings("ConfigFile");
	if (settings.errorOccurred)
		return;
	Array3D_d u(settings.NI, settings.NJ, settings.NK);
	uint counter{0};
	for(uint i{0}; i<settings.NI; ++i)
	{
		for(uint j{0}; j<settings.NJ; ++j)
		{
			for(uint k{0}; k<settings.NK; ++k)
			{
				u(i,j,k) = counter;
				++counter;
			}
		}
	}
	cout << "u(0,0,0) = " << u(0,0,0) << endl;
	cout << "u(0,0,1) = " << u(0,0,1) << endl;
	cout << "u(0,1,0) = " << u(0,1,0) << endl;
	cout << "u(1,0,0) = " << u(1,0,0) << endl;

	cout << "Width: "  << u.getLength()  << endl;
	cout << "Height: " << u.getWidth() << endl;
	cout << "Depth: "  << u.getHeight()  << endl;
	cout << "Size: "   << u.getSize()   << endl;
}

void foo(const int& i)
{
	cout << i << endl;
}

void testPassRvalAsRef()
{
	foo(1+1);
}




#endif /* SRC_TESTS_H_ */
