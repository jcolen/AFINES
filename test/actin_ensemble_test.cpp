#include "actin_ensemble.h"
#include "link_ensemble.h"
#include "Link.h"
#define BOOST_TEST_MODULE actin_ensemble_test
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE( constructors_test )
{
} 

BOOST_AUTO_TEST_CASE( friction_test )
{
}

BOOST_AUTO_TEST_CASE( force_test)
{
}

BOOST_AUTO_TEST_CASE( direction_test)
{
    actin a(0.25, 0.5, 0.75, 1, 2, 2, 4, 4, 0.2);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( ae.get_direction(0)[0] , 0.731689 , tol); 
    BOOST_CHECK_CLOSE( ae.get_direction(0)[1] , 0.681639 , tol); 
}

BOOST_AUTO_TEST_CASE( start_end_test)
{
    actin a(0, 0, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double tol = 0.001; //% 
    BOOST_CHECK_CLOSE( ae.get_start(0)[0], -0.500000, tol); 
    BOOST_CHECK_CLOSE( ae.get_start(0)[1],  0.000000, tol); 
    BOOST_CHECK_CLOSE( ae.get_end(0)[0],  0.500000, tol); 
    BOOST_CHECK_CLOSE( ae.get_end(0)[1],  0.000000, tol); 
    
}

BOOST_AUTO_TEST_CASE( get_intpoint_test)
{
    actin a(0, 0, 0, 1, 0, 0, 0, 0, 0);
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a, 0);
    double x = 0.1, y = 1;
    double tol = 0.001;
    double * intpoint = ae.get_intpoints(0,x,y);
    BOOST_CHECK_CLOSE( intpoint[0], 0.1, tol);
    BOOST_CHECK_CLOSE( intpoint[1], 0, tol);
    double ipx = intpoint[0];
    double ipy = intpoint[1];
    delete[] intpoint;
    BOOST_CHECK_CLOSE( ipx, 0.1, tol);
    BOOST_CHECK_CLOSE( ipy, 0, tol);

}

BOOST_AUTO_TEST_CASE( connect_polymers_test )
{
    actin a1(-1, 0, 0, 1, 0, 0, 0, 0, 0);
    actin a2(1, 0, 0, 1, 0, 0, 0, 0, 0);
    double lnk_len = 1;
    double kl = 100;
    double kb = 2;
    
    link_ensemble * lks = new link_ensemble();
    actin_ensemble ae = actin_ensemble();
    ae.add_monomer(a1, 0);
    ae.add_monomer(a2, 0);
    ae.connect_polymers(lks, lnk_len, kl, kb, "yellow");
    
    BOOST_CHECK_EQUAL(1,1);
    
    //Check that a link exists:
    

}

// EOF
