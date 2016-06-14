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

	jif3D::DATA_STRUCT OutStruct;

	const size_t nsno = 10 + rand() % 100;
	const size_t nrno = 10 + rand() % 100;
	const size_t ntcalc = 10 + rand() % 100;
	const size_t nlshots = 10 + rand() % 100;

	OutStruct.sno.assign(nsno,0.0);
	OutStruct.rno.assign(nrno,0.0);
	OutStruct.tcalc.assign(ntcalc,0.0);
	OutStruct.lshots.assign(nlshots,0.0);
	OutStruct.ndata_seis = rand();
	OutStruct.ndata_seis_act = rand();

	std::generate(OutStruct.sno.begin(),OutStruct.sno.end(),rand);
	std::generate(OutStruct.rno.begin(),OutStruct.rno.end(),rand);
	std::generate(OutStruct.tcalc.begin(),OutStruct.tcalc.end(),drand48);
	std::generate(OutStruct.lshots.begin(),OutStruct.lshots.end(),rand);

	std::ofstream ofs("data_struct.out");
	{
		boost::archive::text_oarchive oa(ofs);
		oa << OutStruct;
	}
	ofs.flush();

	std::ifstream ifs("data_struct.out");
	boost::archive::text_iarchive ia(ifs);
	jif3D::DATA_STRUCT InStruct;
	ia >> InStruct;

	for (size_t i = 0; i < nsno; ++i)
	{
		BOOST_CHECK_EQUAL(OutStruct.sno[i],InStruct.sno[i]);
	}
	for (size_t i = 0; i < nrno; ++i)
	{
		BOOST_CHECK_EQUAL(OutStruct.rno[i],InStruct.rno[i]);
	}
	for (size_t i = 0; i < ntcalc; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.tcalc[i],InStruct.tcalc[i],1e-4);
	}
	for (size_t i = 0; i < nlshots; ++i)
	{
		BOOST_CHECK_EQUAL(OutStruct.lshots[i],InStruct.lshots[i]);
	}
	BOOST_CHECK_EQUAL(OutStruct.ndata_seis,InStruct.ndata_seis);
	BOOST_CHECK_EQUAL(OutStruct.ndata_seis_act,InStruct.ndata_seis_act);
}

BOOST_AUTO_TEST_CASE (rp_struct_test)
{

	jif3D::RP_STRUCT OutStruct;

	const size_t nx = 10 + rand() % 100;
	const size_t ny = 10 + rand() % 100;
	const size_t nz = 10 + rand() % 100;
	const size_t nlen = 10 + rand() % 100;
	const size_t nele = 10 + rand() % 100;

	OutStruct.nray = rand();
	std::generate_n(std::back_inserter(OutStruct.x),nx,drand48);
	std::generate_n(std::back_inserter(OutStruct.y),ny,drand48);
	std::generate_n(std::back_inserter(OutStruct.z),nz,drand48);
	std::generate_n(std::back_inserter(OutStruct.len),nlen,drand48);
	std::generate_n(std::back_inserter(OutStruct.ele),nele,rand);

	std::ofstream ofs("rp_struct.out");
	{
		boost::archive::text_oarchive oa(ofs);
		oa << OutStruct;
	}
	ofs.flush();

	std::ifstream ifs("rp_struct.out");
	boost::archive::text_iarchive ia(ifs);
	jif3D::RP_STRUCT InStruct;
	ia >> InStruct;
	BOOST_CHECK_EQUAL(OutStruct.nray,InStruct.nray);
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
	for (size_t i = 0; i < nlen; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.len[i],InStruct.len[i],1e-4);
	}
	for (size_t i = 0; i < nele; ++i)
	{
		BOOST_CHECK_EQUAL(OutStruct.ele[i],InStruct.ele[i]);
	}

}

BOOST_AUTO_TEST_CASE (ray_result_test)
{
	jif3D::RayResult OutStruct;

	const size_t ntcalc = 10 + rand() % 100;
	const size_t nraypath = 10 + rand() % 100;

	std::generate_n(std::back_inserter(OutStruct.tcalc),ntcalc,drand48);
	for (size_t i = 0; i < nraypath; ++i)
	{
		jif3D::RP_STRUCT ray;
		const size_t nx = 10 + rand() % 100;
		const size_t ny = 10 + rand() % 100;
		const size_t nz = 10 + rand() % 100;
		const size_t nlen = 10 + rand() % 100;
		const size_t nele = 10 + rand() % 100;

		ray.nray = rand();
		std::generate_n(std::back_inserter(ray.x),nx,drand48);
		std::generate_n(std::back_inserter(ray.y),ny,drand48);
		std::generate_n(std::back_inserter(ray.z),nz,drand48);
		std::generate_n(std::back_inserter(ray.len),nlen,drand48);
		std::generate_n(std::back_inserter(ray.ele),nele,rand);

		OutStruct.raypath.push_back(ray);
	}

	std::ofstream ofs("ray_result.out");
	{
		boost::archive::text_oarchive oa(ofs);
		oa << OutStruct;
	}
	ofs.flush();

	std::ifstream ifs("ray_result.out");
	boost::archive::text_iarchive ia(ifs);
	jif3D::RayResult InStruct;
	ia >> InStruct;

	for (size_t i = 0; i < ntcalc; ++i)
	{
		BOOST_CHECK_CLOSE(OutStruct.tcalc[i],InStruct.tcalc[i],1e-4);
	}
	for (size_t i = 0; i < nraypath; ++i)
	{
		BOOST_CHECK_EQUAL(OutStruct.raypath[i].nray,InStruct.raypath[i].nray);
		for (size_t j = 0; j < InStruct.raypath[i].x.size(); ++j)
		{
			BOOST_CHECK_CLOSE(OutStruct.raypath[i].x[j],InStruct.raypath[i].x[j],1e-4);
		}
		for (size_t j = 0; j < InStruct.raypath[i].y.size(); ++j)
		{
			BOOST_CHECK_CLOSE(OutStruct.raypath[i].y[j],InStruct.raypath[i].y[j],1e-4);
		}
		for (size_t j = 0; j < InStruct.raypath[i].z.size(); ++j)
		{
			BOOST_CHECK_CLOSE(OutStruct.raypath[i].z[j],InStruct.raypath[i].z[j],1e-4);
		}
		for (size_t j = 0; j < InStruct.raypath[i].len.size(); ++j)
		{
			BOOST_CHECK_CLOSE(OutStruct.raypath[i].len[j],InStruct.raypath[i].len[j],1e-4);
		}
		for (size_t j = 0; j < InStruct.raypath[i].ele.size(); ++j)
		{
			BOOST_CHECK_EQUAL(OutStruct.raypath[i].ele[j],InStruct.raypath[i].ele[j]);
		}
	}


}

BOOST_AUTO_TEST_SUITE_END()
