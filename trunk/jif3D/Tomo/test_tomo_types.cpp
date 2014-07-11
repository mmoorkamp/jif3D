#define BOOST_TEST_MODULE TomoTypes test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include "tomo_types.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

BOOST_AUTO_TEST_SUITE( Tomo_Types_Test_Suite )

BOOST_AUTO_TEST_CASE (grid_struct_test)
{
	jif3D::GRID_STRUCT OutStruct;

	OutStruct.h = drand48();
	OutStruct.nx = rand();
	OutStruct.ny = rand();
	OutStruct.nz = rand();
	const size_t nslow = 10 + rand() % 100;
	OutStruct.slow.assign(nslow,0.0);
	std::generate(OutStruct.slow.begin(),OutStruct.slow.end(),drand48);

	std::ofstream ofs("grid_struct.out");
	{
		boost::archive::text_oarchive oa(ofs);
		oa << OutStruct;
	}
	ofs.flush();
	std::ifstream ifs("grid_struct.out");
	boost::archive::text_iarchive ia(ifs);
	jif3D::GRID_STRUCT InStruct;
	ia >> InStruct;
	BOOST_CHECK_CLOSE(OutStruct.h,InStruct.h,1e-4);
	BOOST_CHECK_EQUAL(OutStruct.nx,InStruct.nx);
	BOOST_CHECK_EQUAL(OutStruct.ny,InStruct.ny);
	BOOST_CHECK_EQUAL(OutStruct.nz,InStruct.nz);
	for (size_t i = 0; i < nslow; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.slow[i],InStruct.slow[i],1e-4);
	}

}

BOOST_AUTO_TEST_CASE (geometry_test)
{
	jif3D::GEOMETRY OutStruct;

	const size_t nx = 10 + rand() % 100;
	const size_t ny = 10 + rand() % 100;
	const size_t nz = 10 + rand() % 100;
	OutStruct.x.assign(nx,0.0);
	OutStruct.y.assign(ny,0.0);
	OutStruct.z.assign(nz,0.0);
	OutStruct.nrec = rand();
	OutStruct.nshot = rand();

	std::generate(OutStruct.x.begin(),OutStruct.x.end(),drand48);
	std::generate(OutStruct.y.begin(),OutStruct.y.end(),drand48);
	std::generate(OutStruct.z.begin(),OutStruct.z.end(),drand48);

	std::ofstream ofs("geometry.out");
	{
		boost::archive::text_oarchive oa(ofs);
		oa << OutStruct;
	}
	ofs.flush();

	std::ifstream ifs("geometry.out");
	boost::archive::text_iarchive ia(ifs);
	jif3D::GEOMETRY InStruct;
	ia >> InStruct;

	for (size_t i = 0; i < nx; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.x[i],InStruct.x[i],1e-4);
	}
	for (size_t i = 0; i < ny; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.y[i],InStruct.y[i],1e-4);
	}
	for (size_t i = 0; i < nz; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.z[i],InStruct.z[i],1e-4);
	}
	BOOST_CHECK_EQUAL(OutStruct.nrec,InStruct.nrec);
	BOOST_CHECK_EQUAL(OutStruct.nshot,InStruct.nshot);
}

BOOST_AUTO_TEST_CASE (data_struct_test)
{

}

BOOST_AUTO_TEST_CASE (rp_struct_test)
{

}

BOOST_AUTO_TEST_CASE (ray_result_test)
{

}

BOOST_AUTO_TEST_SUITE_END()
