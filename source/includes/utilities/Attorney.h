//
//  Attorney.h
//  CSMP_ReservoirSimulator
//

#ifndef ATTORNEY_H
#define ATTORNEY_H

namespace csmp {

/**

@brief Attorney design pattern gives access to protected / private member
variables and methods.


@section motivation Motivation

Unit tests may require access to internal variables and methods to test the
basic internal works of a class.


@section design Design Intent

This version is templated to
     1. avoid friend declaration in the class being tested. Cluttering a class
     with friend declarations which are not relevant to the way the class is
     used in the code is not helpful to the user.
     2. give selective access to members and methods. A friend declaration gives
     access to all protected and private members and methods. This templated
     version requires explicit declaration of which members become accessible.


@section example Application Example

Declare an Attorney of the class to be tested and declare members and methods to
be made accessible like so:

@code
template<uint32_t dim>
class PDE_Integrator_Attorney : public Attorney<class PDE_Integrator<dim>> {
public:
    using Attorney<PDE_Integrator<dim>>::Attorney; // Inherit Attorney constructor
    using PDE_Integrator<dim>::EstablishMatrixSetup;
    using PDE_Integrator<dim>::lhs_operators_;
};
@endcode

Then, in your test, create an Attorney and use it like so:

@code
PDE_Integrator<2U,Element> pde_integrator( ... );
PDE_Integrator_Attorney<2U> attorney( pde_integrator );
LHS_FixedValueMatrix<2U>  lhs( ... );
attorney.Add(&lhs);
attorney.EstablishMatrixSetup( ... ); // <- Accessing protected method
_test( attorney.lhs_operators_.size() == 1 ); // <- Accessing protected member variable
@endcode

*/
template<class Client>
class Attorney : public Client {
public:
    Attorney( Client& client ) : client_( client ) {
    }
private:
    Client& client_;
};
} // end csmp

#endif /* defined(ATTORNEY_H) */

