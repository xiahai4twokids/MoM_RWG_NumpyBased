#include <iostream>
#include <boost/array.hpp>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "boost/lambda/lambda.hpp"
#include "boost/function.hpp"
#include "boost/python.hpp"
#include "boost/python/numpy/ndarray.hpp"

using namespace boost::numeric::ublas;
using namespace boost;
using namespace boost::lambda;
namespace p = boost::python;
namespace np = boost::python::numpy;
typedef std::complex<float> complex;




struct Result_RankRevealingMethod_tradition
{
	shared_ptr<
		matrix<complex>
	> L;
	shared_ptr<
		matrix<complex>
	> U;
	float res;
	int k;

	Result_RankRevealingMethod_tradition(shared_ptr<matrix<complex>> LL, shared_ptr<matrix<complex>> UU, float res, int k):
	L(LL), U(UU), res(res), k(k){}
	Result_RankRevealingMethod_tradition(){}
};

class matrix_defined
{	
	shared_ptr<matrix<float>> field;
	shared_ptr<matrix<float>> source;
	float wavenumber;	
public:
	matrix_defined()
	{
		this->field = shared_ptr<matrix<float>>(
			new matrix<float>(4,3));
		this->source = shared_ptr<matrix<float>>(
			new matrix<float>(4,3));
		for (int i=0; i<this->field->size1(); i++)
		{
			this->field->insert_element(i,0,i);
			this->field->insert_element(i,1,0);
			this->field->insert_element(i,2,0);
		}
		for (int i=0; i<this->source->size1(); i++)
		{
			this->source->insert_element(i,0,i+10);
			this->source->insert_element(i,1,0);
			this->source->insert_element(i,2,0);
		}
		this->wavenumber = float(2.*3.14159);
	}

	matrix_defined(shared_ptr<matrix<float>> field, shared_ptr<matrix<float>> source, float wavenumber):field(field),source(source),wavenumber(wavenumber){}
	shared_ptr<matrix<complex>> print()
	{
		shared_ptr<matrix<complex>> t;
	       	shared_ptr<vector<complex>> temp;
		try{ 
			t =  shared_ptr<matrix<complex>>(
				new matrix<complex>(this->field->size1(), this->source->size1()) 
				); 
		}
	       	catch(...) { throw "print 1"; }
		for (int i=0; i<t->size1(); i++)
		{
			try{ temp = this->getRow(i+1); }
			catch(...) { throw "print 2"; } 
			for(int j=0; j<t->size2(); j++)
			try{ (*t)(i,j) = (*temp)(j); }
			catch(...) { throw "print 3"; }
		}
		return t;
	}
    complex kernel(float r){	return exp(complex(0,-this->wavenumber*r));  	}

	int size1()
	{
		return this->field->size1();
	}
	int size2()
	{
		return this->source->size1();
	}


	shared_ptr<vector<complex>> getCol(int k)
	{ 
		shared_ptr<vector<complex>> a(
			new vector<complex>(this->field->size1())
			);
		vector<complex> source (
			row( *(this->source), k-1	)
			);
		for (int i=0; i<a->size(); i++)
		{
			vector<complex> field ( 
					row(*(this->field), i)
				        );
			field -= source;
			float dist = norm_2(field);
			a->insert_element(i, this->kernel(dist));
		}
		// std::cout << a->size() << std::endl;
		return a; 
	}
	shared_ptr<vector<complex>> getRow(int k)
	{ 
		shared_ptr<vector<complex>> a(
			new vector<complex>(this->source->size1())
			);
		vector<complex> field (
			       	row( *(this->field), k-1)
				);
		for (int i=0; i<a->size(); i++)
		{			
			vector<complex> source (
				       	row(*(this->source), i)
					);
			source -= field;
			float dist = norm_2(source);
			a->insert_element(i, this->kernel(dist));
		}
		return a; 
	}
};

class RankRevealingMethod_tradition
{
	float threshold_NULL;
public:
	RankRevealingMethod_tradition(){ this->threshold_NULL = 0.00001f;}
public:
	Result_RankRevealingMethod_tradition LUDec (
			matrix_defined _matrix, int rows, int columns, 
			float threshold_remain, int threshold_rank)
	{
		matrix<complex> B(zero_matrix<complex>(rows, columns));

		int rank = 0;
		if(-1 == threshold_rank ){ rank = columns; }
		else{ rank = columns; }
		matrix<complex> L(zero_matrix<complex>(columns, rows));
		matrix<complex> U(zero_matrix<complex>(rows, columns));
		vector<complex> l(rows);
		vector<complex> u(columns);

		float Z_norm = 0.;		
		int j = 0;
		float max = 0.;
		complex max_lj(0,0);
		float addMatrix_norm = 0;
		float res = 100;

	        int k=1;
		for ( k=1; k<columns+1; k++)
		{
			try{  l = *(_matrix.getCol(k)) - column(B,k-1);	}
			catch(...){ 
				std::cout << "LUDec0 0\n";
			       	throw "Exception 0";  
			}			
			
			j = 0; // j = np.argmax(np.abs(l))+1
			max = 0.;
			max_lj = 0.;
			for (int i=0; i<l.size();i++)
			{
				try{		
					if (max<abs(l(i))){j=i;max=abs(l(i));}						
				}catch(...){
				        std::cout << "LUDec0 1\n";
					throw "Exception 1";
				}
			}
			j += 1;
			if (max < this->threshold_NULL) 
			{/*std::cout << "max="<<max<<"break\n"; */k--; break; }
			try{ max_lj = l(j-1);	}
			catch(...){
				std::cout << "LUDEC0 2\n";
				throw "Exception 2";	
			}
			if ( abs(max_lj)<threshold_NULL )
			{/*std::cout << "max:" << max << "\tcontinue\n";*/  continue;	}
			
			try{  u = *(_matrix.getRow(j)) - row(B,j-1); }
			catch(...){
				std::cout << "LUDEC0 3\n";
				throw "Exception 3";	
			}
			l = l/max_lj;

			matrix<complex> l_mat(l.size(),1);
			matrix<complex> u_mat(1,u.size());
			try{
				for (int i=0;i<l.size();i++) l_mat(i,0)=l(i);
				for (int i=0;i<u.size();i++) u_mat(0,i)=u(i);
			}catch(...){
				std::cout << "LUDEC0 4\n";	
				throw "Exception 4";	
			}
			try{ B = B+prod(l_mat,u_mat); }
			catch(...){	
				std::cout << "LUDEC0 5\n";
				throw "Exception 5";
			}

			addMatrix_norm = norm_2(l)*norm_2(u);
			Z_norm = Z_norm + addMatrix_norm*addMatrix_norm;
			for (int ii=0; ii<k-1;ii++)	{
				try{
					Z_norm = Z_norm +2.*abs(
						inner_prod( row (L, ii),l )*inner_prod( row (U, ii),u )
						);
				}catch(...){	
					std::cout << "LUdec0 6\n";
					throw "Exception 6";	
				}
			}
			
			res = addMatrix_norm/sqrt(Z_norm);
			for (int i=0; i<l.size(); i++)
				L.insert_element (k-1,i,l(i));
			for (int i=0; i<u.size(); i++)
				U.insert_element(k-1,i,u(i));
		        
			if (res<=threshold_remain || k>=rank)
			{/*std::cout << "k="<<k << std::endl;*/ break; }			
		}

		shared_ptr<	matrix<complex>> LL( 
			new matrix<complex>(
			zero_matrix<complex>(L.size2(),k)
			)
			);
		shared_ptr<	matrix<complex>> UU( 
			new matrix<complex>(
			zero_matrix<complex>(k,U.size2())
			)
			);
		for(int i=0; i<LL->size1(); i++)
			for (int j=0; j<LL->size2(); j++)
				try{ (*LL)(i,j) = L(j,i);}catch(...){std::cout << "LUDEC0 7\n";throw "Exception 7"; }
		for(int i=0; i<UU->size1(); i++)
			for (int j=0; j<UU->size2(); j++)
				try{(*UU)(i,j) = U(i,j);}catch(...){
					std::cout << i<< "  " << j << std::endl;
					std::cout << UU->size1() << "  " << UU->size2() << std::endl;
					std::cout << U.size1() << "  " << U.size2() << std::endl;
					std::cout << "LUDEC0 8\n"; throw "Exception 8";
			       	}
		
		return Result_RankRevealingMethod_tradition(
			LL,UU,res,k
			);
	}

};

///////////////////////////////////////////////////////////////
p::list matrixC2P_complex(shared_ptr<matrix<complex>> mat)
{
	p::list tempMat;
	for (int i=0; i<mat->size1(); i++)
	{
		p::list row;
		for (int j=0; j<mat->size2();j++)
		{
			try{
				row.append( (*mat)(i,j));
			}catch(...)
			{ 
				std::cout << "matrixC2P\n";
				throw "matrixC2P"; 
			}
		}
		tempMat.append(row);
	}
	return tempMat;
}

shared_ptr<matrix<float>> matrixP2C_float(p::list mat )
{
	shared_ptr<matrix<float>> tempMat( new 
			matrix<float>(p::len(mat), p::len(mat[0])));
	shared_ptr<vector<float>> tempVec;
	for (int i=0; i<tempMat->size1(); i++)
		for (int j=0; j<tempMat->size2(); j++)
		{
			try{
				(*tempMat)(i,j) = p::extract<float>(mat[i][j]);
			}catch(...)
			{
				std::cout << "matrixP2C\n";
				throw "matrixP2C";
			}
		}
	return tempMat;
}
class RankRevealingMethod_tradition_python
{
public:
  	RankRevealingMethod_tradition_python(){}
	p::list LUDec(p::list field, p::list source, p::object wavenumber,
			p::object rows, p::object columns, 
			p::object threshold_remain, p::object threshold_rank
			)
	{
		Result_RankRevealingMethod_tradition result_C;
		try{ 
			result_C = RankRevealingMethod_tradition().LUDec(
				matrix_defined(matrixP2C_float(field),matrixP2C_float(source),p::extract<float>(wavenumber)),
				p::extract<int>(rows), p::extract<int>(columns),
				p::extract<float>(threshold_remain), p::extract<float>(threshold_rank)
				);
		}catch(...)
		{ 
			std::cout << "LUDEc0\n";
			throw "RRM_tra_python";
		}
		p::list result_P;
		try{
			result_P.append(matrixC2P_complex(result_C.L));
			result_P.append(matrixC2P_complex(result_C.U));
			result_P.append(result_C.res);
			result_P.append(result_C.k);
		}catch(...)
		{
			std::cout << "LUDec1\n";
			throw "RRM_tra_python 1";
		}

		return result_P; 

	}
};
