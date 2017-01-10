#include "GradientGrid.h"
#include <list>
#include "../constants.h"

void GradientGrid::clear()
{
	dzdxOut.clear(); 
	dzdyOut.clear(); 
	magOut.clear(); 
	slopeOut.clear(); 
	aspectOut.clear();
	existFlag = false;
}

void GradientGrid::calc_GradientGrid(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn)
{
	//*******************************************************************************************
	//A. Prepare inputs x, y, z
	//*******************************************************************************************
	std::list<double> xIn2;
	xIn2.assign(xIn.begin(),xIn.end());
	std::list<double> yIn2;
	yIn2.assign(yIn.begin(),yIn.end());

	xIn2.sort();
	yIn2.sort();
	xIn2.unique();
	yIn2.unique();

	int nc = (const int)xIn2.size();
	int nr = (const int)yIn2.size();
	int nc2 = nc+2;
	int nr2 = nr+2;
	
	dgrid x = dgrid(xIn);
	dgrid y = dgrid(yIn);
	dgrid z = dgrid(zIn);
	x.resize(nr,nc);
	y.resize(nr,nc);	
	z.resize(nr,nc);

	//*******************************************************************************************
	//B. Pre-initialize storage arrays
	//*******************************************************************************************
	//Gradient attributes
	dzdxOut = dgrid(nr,nc,0.00);
	dzdyOut = dgrid(nr,nc,0.00);
	slopeOut = dgrid(nr,nc,0.00);
	aspectOut = dgrid(nr,nc,0.00);
	magOut = dgrid(nr,nc,0.00);

	//Temps
	dgrid z_V2(nr2,nc2,0.00);
	dgrid x_V2(nr2,nc2,0.00);
	dgrid y_V2(nr2,nc2,0.00);
	dgrid dzdxOut2(nr2,nc2,0.00);
	dgrid dzdyOut2(nr2,nc2,0.00);
	dgrid dxdrow(nr2,nc2,0.00);
	dgrid dydrow(nr2,nc2,0.00);
	dgrid dxdcol(nr2,nc2,0.00);
	dgrid dydcol(nr2,nc2,0.00);
	dgrid div(nr2,nc2,0.00);
	dgrid ax(nr2,nc2,0.00);
	dgrid ay(nr2,nc2,0.00);
	dgrid a(nr2,nc2,0.00);
	dgrid bx(nr2,nc2,0.00);
	dgrid by(nr2,nc2,0.00);
	dgrid b(nr2,nc2,0.00);
	dgrid cx(nr2,nc2,0.00);
	dgrid cy(nr2,nc2,0.00);
	dgrid c(nr2,nc2,0.00);
	dgrid ex(nr2,nc2,0.00);
	dgrid ey(nr2,nc2,0.00);
	dgrid e(nr2,nc2,0.00);
	dgrid gx(nr2,nc2,0.00);
	dgrid gy(nr2,nc2,0.00);
	dgrid g(nr2,nc2,0.00);
	dgrid lx(nr2,nc2,0.00);
	dgrid ly(nr2,nc2,0.00);
	dgrid l(nr2,nc2,0.00);
	dgrid mx(nr2,nc2,0.00);
	dgrid my(nr2,nc2,0.00);
	dgrid m(nr2,nc2,0.00);
	dgrid nx(nr2,nc2,0.00);
	dgrid ny(nr2,nc2,0.00);
	dgrid n(nr2,nc2,0.00);
	dgrid dzdrow(nr2,nc2,0.00);
	dgrid dzdcol(nr2,nc2,0.00);
	dgrid colang(nr2,nc2,0.00);
	dgrid rowang(nr2,nc2,0.00);
	dgrid rowdist(nr2,nc2,0.00);
	dgrid coldist(nr2,nc2,0.00);
	dgrid dzdx(nr2,nc2,0.00);
	dgrid dzdy(nr2,nc2,0.00);
	nc = (const int)xIn2.size()-1;
	nr = (const int)yIn2.size()-1;
	nr2--;
	nc2--;
		
#pragma region -- Gradient Calculations

	int ii, jj;
	for(int i = 0; i <= nr; i++)
	{
		for(int j=0;j<=nc;j++)
		{
			ii = i+1;
			jj = j+1;
			
			//%% I. Assign data to "V2" arrays.
			//%A. Populate center portions
			x_V2(ii,jj) = x(i,j);
			y_V2(ii,jj) = y(i,j);
			z_V2(ii,jj) = z(i,j);
		}
	}

	for(int j=0;j<=nc;j++)
	{	
		//%B. Get edges of V2 matrices
		x_V2(0,j+1) = 2*x(0,j)-x(1,j);
		x_V2(nr2,j+1) = 2*x(nr,j)-x(nr-1,j);
		y_V2(0,j+1) = 2*y(0,j)-y(1,j);
		y_V2(nr2,j+1) = 2*y(nr,j)-y(nr-1,j);
		z_V2(0,j+1) = 2*z(0,j)-z(1,j);
		z_V2(nr2,j+1) = 2*z(nr,j)-z(nr-1,j);
	}
	
	for(int i=0;i<=nr;i++)
	{
		//%B. Get edges of V2 matrices
		x_V2(i+1,0)=2*x(i,0)-x(i,1);
		x_V2(i+1,nc2) = 2*x(i,nc)-x(i,nc-1);
		y_V2(i+1,0)=2*y(i,0)-y(i,1);
		y_V2(i+1,nc2) = 2*y(i,nc)-y(i,nc-1);
		z_V2(i+1,0)=2*z(i,0)-z(i,1);
		z_V2(i+1,nc2) = 2*z(i,nc)-z(i,nc-1);
	}

	//%C. Get Corners
	x_V2(0,0) = (2*x_V2(1,0)-x_V2(2,0))/2 + (2*x_V2(0,1)-x_V2(0,2))/2;
	x_V2(0,nc2) = (2*x_V2(0,nc2-1)-x_V2(0,nc2-2))/2 + (2*x_V2(1,nc2)-x_V2(2,nc2))/2;
	x_V2(nr2,0) = (2*x_V2(nr2-1,0)-x_V2(nr2-2,0))/2 + (2*x_V2(nr2,1)-x_V2(nr2,2))/2;
	x_V2(nr2,nc2) = (2*x_V2(nr2-1,nc2)-x_V2(nr2-2,nc2))/2 + (2*x_V2(nr2,nc2-1)-x_V2(nr2,nc2-2))/2;
	
	y_V2(0,0) = (2*y_V2(1,0)-y_V2(2,0))/2 + (2*y_V2(0,1)-y_V2(0,2))/2;
	y_V2(0,nc2) = (2*y_V2(0,nc2-1)-y_V2(0,nc2-2))/2 + (2*y_V2(1,nc2)-y_V2(2,nc2))/2;
	y_V2(nr2,0) = (2*y_V2(nr2-1,0)-y_V2(nr2-2,0))/2 + (2*y_V2(nr2,1)-y_V2(nr2,2))/2;
	y_V2(nr2,nc2) = (2*y_V2(nr2-1,nc2)-y_V2(nr2-2,nc2))/2 + (2*y_V2(nr2,nc2-1)-y_V2(nr2,nc2-2))/2;

	z_V2(0,0) = (2*z_V2(1,0)-z_V2(2,0))/2 + (2*z_V2(0,1)-z_V2(0,2))/2;
	z_V2(0,nc2) = (2*z_V2(0,nc2-1)-z_V2(0,nc2-2))/2 + (2*z_V2(1,nc2)-z_V2(2,nc2))/2;
	z_V2(nr2,0) = (2*z_V2(nr2-1,0)-z_V2(nr2-2,0))/2 + (2*z_V2(nr2,1)-z_V2(nr2,2))/2;
	z_V2(nr2,nc2) = (2*z_V2(nr2-1,nc2)-z_V2(nr2-2,nc2))/2 + (2*z_V2(nr2,nc2-1)-z_V2(nr2,nc2-2))/2;

	//%% II. Get data points for third order technique
	for (int i=0;i<=nr2;i++)//nr
	{
		for (int j=0;j<=nc2;j++)
		{
		    div(i,j) = 0.0;
			//%A. For a
			if ((i-1) < 0 || (j-1) < 0)
			{
				ax(i,j) = 0.0;
				ay(i,j) = 0.0;
				a(i,j) = 0.0;
			}
			else
			{
				ax(i,j) = x_V2(i-1,j-1);
				ay(i,j) = y_V2(i-1,j-1);
				a(i,j) = z_V2(i-1,j-1);
				div(i,j)=div(i,j)+1;
			}
			//B. For b
			if (i-1 < 0)
			{
				bx(i,j) = 0.0;
				by(i,j) = 0.0;
				b(i,j) = 0.0;
			}
			else
			{
				bx(i,j) = x_V2(i-1,j);
				by(i,j) = y_V2(i-1,j);
				b(i,j) = z_V2(i-1,j);
				div(i,j)=div(i,j)+1;
			}
			//%C. For c
			if ((i-1) < 0 || (j+1) == nc2+1)
			{
				cx(i,j) = 0.0;
				cy(i,j) = 0.0;
				c(i,j) = 0.0;
			}
			else
			{
				cx(i,j) = x_V2(i-1,j+1);
				cy(i,j) = y_V2(i-1,j+1);
				c(i,j) = z_V2(i-1,j+1);
				div(i,j)=div(i,j)+1;
			}
			//%D. For e
			if ((j-1) < 0)
			{
				ex(i,j) = 0.0;
				ey(i,j) = 0.0;
				e(i,j) = 0.0;
			}
			else
			{
				ex(i,j) = x_V2(i,j-1);
				ey(i,j) = y_V2(i,j-1);
				e(i,j) = z_V2(i,j-1);
				div(i,j)=div(i,j)+1;
			}
			//%E. For g
			if ((j+1) == nc2+1)
			{
				gx(i,j) = 0.0;
				gy(i,j) = 0.0;
				g(i,j) = 0.0;
			}
			else
			{
				gx(i,j) = x_V2(i,j+1);
				gy(i,j) = y_V2(i,j+1);
				g(i,j) = z_V2(i,j+1);
				div(i,j)=div(i,j)+1;
			}
			//%F. For l
			if ((i+1) == nr2+1 || (j-1) < 0)
			{
				lx(i,j) = 0.0;
				ly(i,j) = 0.0;
				l(i,j) = 0.0;
			}
			else
			{
				lx(i,j) = x_V2(i+1,j-1);
				ly(i,j) = y_V2(i+1,j-1);
				l(i,j) = z_V2(i+1,j-1);
				div(i,j)=div(i,j)+1;
			}
			//%G. For m
			if ((i+1) == nr2+1)
			{
				mx(i,j) = 0.0;
				my(i,j) = 0.0;
				m(i,j) = 0.0;
			}
			else
			{
				mx(i,j) = x_V2(i+1,j);
				my(i,j) = y_V2(i+1,j);
				m(i,j) = z_V2(i+1,j);
				div(i,j)=div(i,j)+1;
			}
			//%H. For n
			if ((i+1) == nr2+1 || (j+1) == nc2+1)
			{
				nx(i,j) = 0.0;
				ny(i,j) = 0.0;
				n(i,j) = 0.0;
			}
			else
			{
				nx(i,j) = x_V2(i+1,j+1);
				ny(i,j) = y_V2(i+1,j+1);
				n(i,j) = z_V2(i+1,j+1);
				div(i,j)=div(i,j)+1;
			}

			//%Roman I. Calculate east-west and north-south gradients 
			//%(size is 1 under normal conditions)
			dydrow(i,j) = ((cy(i,j) + 2 * gy(i,j) + ny(i,j)) - (ay(i,j)+ 2 * ey(i,j) + ly(i,j))) / double(div(i,j));
			dydcol(i,j) = ((ay(i,j) + 2 * by(i,j) + cy(i,j)) - (ly(i,j) + 2 * my(i,j) + ny(i,j))) / double(div(i,j));
			dxdrow(i,j) = ((cx(i,j) + 2 * gx(i,j) + nx(i,j)) - (ax(i,j) + 2 * ex(i,j) + lx(i,j))) / double(div(i,j));
			dxdcol(i,j) = ((ax(i,j) + 2 * bx(i,j) + cx(i,j)) - (lx(i,j) + 2 * mx(i,j) + nx(i,j))) / double(div(i,j));
			dzdrow(i,j) = ((c(i,j) + 2 * g(i,j) + n(i,j)) - (a(i,j) + 2 * e(i,j) + l(i,j))) / double(div(i,j));
			dzdcol(i,j) = ((a(i,j) + 2 * b(i,j) + c(i,j)) - (l(i,j) + 2 * m(i,j) + n(i,j))) / double(div(i,j));
           
			colang(i,j) = atan2(dydcol(i,j),dxdcol(i,j));
			rowang(i,j) = atan2(dydrow(i,j),dxdrow(i,j));
           
			rowdist(i,j) = sqrt(pow(dxdrow(i,j),2) + pow(dydrow(i,j),2));
			coldist(i,j) = sqrt(pow(dxdcol(i,j),2) + pow(dydcol(i,j),2));
           
			dzdx(i,j) =  dzdcol(i,j) / coldist(i,j) * cos(colang(i,j)) + dzdrow(i,j) / rowdist(i,j) * cos(rowang(i,j));
			dzdy(i,j) =  dzdrow(i,j) / rowdist(i,j) * sin(rowang(i,j)) + dzdcol(i,j) / coldist(i,j) * sin(colang(i,j));
			
			//%J. to prevent division by zero
			if (dzdx(i,j) == 0.0)
				dzdx(i,j) = 0.0000001;
			
			//%K. calculate the slope and return
			dzdxOut2(i,j) = dzdx(i,j);
			dzdyOut2(i,j) = dzdy(i,j);
			
			//%% III. Assign data to output arrays
			if((j<=nc+1 && i<=nr+1)&& (j!=0 && i!=0))
			{
				//%A. Gradient in the E-W direction
				dzdxOut(i-1,j-1) = dzdxOut2(i,j);

				//%B. Gradient in the N-S direction
				dzdyOut(i-1,j-1) = dzdyOut2(i,j);

				//%C. Magnitude of the gradient 
				magOut(i-1,j-1) = sqrt(pow(dzdxOut(i-1,j-1),2) + pow(dzdyOut(i-1,j-1),2));
			
				//%D. Upward slope in degrees relevant to a flat bottom
				slopeOut(i-1,j-1) = 180.0/PI*atan(magOut(i-1,j-1));

				//%E. Clockwise direction of steepest ascent, with zero degrees being north.
				aspectOut(i-1,j-1) = 270.0 + 180.0/PI*atan(dzdyOut(i-1,j-1)/(dzdxOut(i-1,j-1))) - 90.0*(dzdxOut(i-1,j-1)/abs(dzdxOut(i-1,j-1)));
			}
		}
	}
	existFlag = true; //Remember we calculated gradients
	#pragma endregion Gradient Calculations
	
	xIn2.clear();
	yIn2.clear();
	x.clear();
	y.clear();
	z.clear();
	z_V2.clear();
	x_V2.clear();
	y_V2.clear();
	dzdxOut2.clear();
	dzdyOut2.clear();
	dxdrow.clear();
	dydrow.clear();
	dxdcol.clear();
	dydcol.clear();
	div.clear();
	ax.clear();
	ay.clear();
	a.clear();
	bx.clear();
	by.clear();
	b.clear();
	cx.clear();
	cy.clear();
	c.clear();
	ex.clear();
	ey.clear();
	e.clear();
	gx.clear();
	gy.clear();
	g.clear();
	lx.clear();
	ly.clear();
	l.clear();
	mx.clear();
	my.clear();
	m.clear();
	nx.clear();
	ny.clear();
	n.clear();
	dzdrow.clear();
	dzdcol.clear();
	colang.clear();
	rowang.clear();
	rowdist.clear();
	coldist.clear();
	dzdx.clear();
	dzdy.clear();
}


