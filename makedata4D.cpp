//make4Ddata.cpp

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <string>


using namespace std;

double dist2(double dx, double dy, double dz, double dw, int xx, int yy, int zz, int ww){
	return (dx-1.0*xx)*(dx-1.0*xx)+(dy-1.0*yy)*(dy-1.0*yy)+(dz-1.0*zz)*(dz-1.0*zz)+(dw-1.0*ww)*(dw-1.0*ww);
}

void print_usage_and_exit(int exit_code) {
	cerr << "Usage: "
	     << "MD4D "
	     << "[options] " << endl
	     << endl
	     << "Options:" << endl
	     << endl
	     << "  --help           print this screen" << endl
	     << "  --nsize <n>      real data size, <n>^4" << endl
	     << "  --ssize <s>      pattern size, once we make a random data of <s>^4, and spread it into <s>^4" << endl
	     << "  --radius <r>     when we spread the pattern into real data, <r> is the blur radius." << endl
	     << " output filename is MD4D_<n>_<s>_<r>.complex"
	     << endl;

	exit(exit_code);
}


int main(int argc, char** argv){

	int n=20;
	int s=10;
	int area_radius = 3;

	for (int i = 1; i < argc; ++i) {
		const string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--nsize") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			n = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--ssize") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			s = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--radius") {
			string parameter = string(argv[++i]);
			size_t next_pos;
			area_radius = stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} 
	}


	random_device rd;	
	mt19937 mt(rd());
	uniform_real_distribution<double> score(0.0,1.0);


	double  data[s][s][s][s];
	for(int x=0; x<s; ++x)for(int y=0; y<s; ++y)for(int z=0; z<s; ++z)for(int w=0; w<s; ++w){
		data[x][y][z][w] = score(mt);
	}


	string n_str=to_string(n);
	string s_str=to_string(s);

	double rate = 1.0*s /n;

	string radius_str = to_string(area_radius);
	double area = 0.5*area_radius;



	string outname = "MD4D_"+n_str+"_"+s_str+"_"+radius_str+".complex";
	ofstream writing_file;
	writing_file.open(outname, ios::out | ios::binary);

	if(!writing_file.is_open()){
		cout << " error: open file for output failed! " << endl;
	}

	int64_t mn = 8067171840;
	writing_file.write((char *) &mn, sizeof( int64_t )); // magic number
	
	int64_t type = 1;
	writing_file.write((char *) &type, sizeof( int64_t )); // type number of PERSISTENCE_DIAGRAM
	
	int64_t totalnum = n*n*n*n;
	writing_file.write((char *) &totalnum, sizeof( int64_t )); // total size of data

	int64_t dim = 4;
	writing_file.write((char *) &dim, sizeof( int64_t )); // dimension

	int64_t numx=n;
	writing_file.write((char *) &numx, sizeof( int64_t )); // size of x
	writing_file.write((char *) &numx, sizeof( int64_t )); // size of y
	writing_file.write((char *) &numx, sizeof( int64_t )); // size of z
	writing_file.write((char *) &numx, sizeof( int64_t )); // size of w
	
	

	for(int x=0; x<n; x++)for(int y=0; y<n; y++)for(int z=0; z<n; z++)for(int w=0; w<n; w++){
		double dx = 0.1 + rate*x;
		double dy = 0.1 + rate*y;
		double dz = 0.1 + rate*z;
		double dw = 0.1 + rate*w;
		int startx = max((int)floor(dx-area),0);
		int starty = max((int)floor(dy-area),0);
		int startz = max((int)floor(dz-area),0);
		int startw = max((int)floor(dw-area),0);
		int endx = min((int)floor(dx+area),s-1);
		int endy = min((int)floor(dy+area),s-1);
		int endz = min((int)floor(dz+area),s-1);
		int endw = min((int)floor(dw+area),s-1);
		double sum1=0, sum2=0;

		for(int xx=startx; xx<=endx; xx++)for(int yy=starty; yy<=endy; yy++)for(int zz=startz; zz<=endz; zz++)for(int ww=startw; ww<=endw; ww++){
			double d2 = dist2(dx,dy,dz,dw,xx,yy,zz,ww);
			sum1 += (1/(1+d2*d2));
			sum2 += (data[xx][yy][zz][ww]/(1+d2*d2));
		}
		double v = sum2 / sum1;
		//make
		//cout << v << " ";
		writing_file.write((char *) &v, sizeof( double )); // birth
	}

	writing_file.close();

	cout << outname << endl;
	return 1;

}