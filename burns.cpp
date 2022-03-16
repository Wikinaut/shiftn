// implementation of line detection as described in http://dx.doi.org/10.1109/TPAMI.1986.4767808 
// Marcus Hebel

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include <iostream>
#include "mystack.h"

using namespace std;
using namespace cv;

const int layer_c_intern = 11;
const int layer_l_intern = 2;
const double mypi = 3.14159265358979323846264338328;

// line extraction parameter
int numpart = 10;	// number of partitions of 360 degrees
int minsize = 46;	// min region size for line (approx. double line lenghth)
int minlen = 23;	// min line length
int mincon = 3;	// minimum contrast of edges
double mingrad = 1.5;	// minimum gradient of edges
double minsup = 0.6;	// minimum portion of pixels supporting the line
int rgb[] = { 30,59,11 }; // rgb portions to convert to gray-scale

Mat theimage[layer_c_intern];	// data structure to store the image channels (r, g, b, gray, low-pass, gradient, lines, tmp1, tmp2, tmp3, tmp4) 
Mat l_layer[layer_l_intern];	// data structure for calculations
long *structure;	// data structure to store the list of detected lines
long structure_set = 0;
int *grayvalues;
int memory_set = 0;
int result_set = 0;
int headersize = 15;
int x_global = 0;
int y_global = 0;

int uc_int(uchar c) {
	int ret = (int)c;
	return ret;
}

uchar i_uchar(int i) {
	uchar c;
	if (i > 255) i = 255;
	if (i < 0) i = 0;
	c = (uchar)i;
	return c;
}

int read_theimage(int id, int &x, int &y, string mfilename) {	// read image using OpenCV
	Mat cvimg;
	int i;
	int j;
	int data;
	int u = 0;
	int v = 0;
	int ret = 0;
	uchar c;
	cvimg = imread(mfilename);
	u = cvimg.cols;
	v = cvimg.rows;
	if ((u == 0) || (v == 0) || (u > 65535) || (v > 65535)) ret = 1;
	if (!ret) {
		for (i = 0; i < layer_c_intern; i++) {
			theimage[i] = Mat::zeros(v, u, CV_8U);
		}
		for (i = 0; i < layer_l_intern; i++) {
			l_layer[i] = Mat::zeros(v, u, CV_32S);
		}
	}
	if (!ret) {
		memory_set = 1;
		for (j = 0; j < v; j++) {
			for (i = 0; i < u; i++) {
				data = 0;
				c = cvimg.at<Vec3b>(j, i)[0];
				theimage[id].at<uchar>(j, i) = c;
				data += rgb[2] * uc_int(c);
				c = cvimg.at<Vec3b>(j, i)[1];
				theimage[id + 1].at<uchar>(j, i) = c;
				data += rgb[1] * uc_int(c);
				c = cvimg.at<Vec3b>(j, i)[2];
				theimage[id + 2].at<uchar>(j, i) = c;
				data += rgb[0] * uc_int(c);
				data /= (rgb[0] + rgb[1] + rgb[2]);
				theimage[id + 3].at<uchar>(j, i) = i_uchar(data);
			}
		}
	}
	x = u;
	y = v;
	cvimg.release();
	return ret;
}

void gauss(int id_old, int id_new) { // low-pass filter of grayscale image
	int i;
	int j;
	int k;
	int l;
	int thesum;
	int thecount;
	int themultiplicator;
	for (j = 0; j < y_global; j++) {
		for (i = 0; i < x_global; i++) {
			thesum = 0;
			thecount = 0;
			for (k = -1; k <= 1; k++) {
				for (l = -1; l <= 1; l++) {
					if ((i + k >= 0) && (i + k < x_global) && (j + l >= 0) && (j + l < y_global)) {
						themultiplicator = 3 - k * k - l * l;
						if ((k == 0) && (l == 0)) themultiplicator += 1;
						thesum += themultiplicator * uc_int(theimage[id_old].at<uchar>(j + l, i + k));
						thecount += themultiplicator;
					}
				}
			}
			theimage[id_new].at<uchar>(j, i) = i_uchar((int)(thesum / thecount));
		}
	}
}

int partitionvalue(double w, double wsteph) {
	int i;
	int ret = 0;
	double wdown;
	wdown = -mypi;
	i = 1;
	while (i <= numpart) {
		if ((w >= wdown) && (w < (wdown + wsteph))) {
			ret = i;
			i = numpart;
		}
		wdown += wsteph;
		i += 1;
	}
	return ret;
}

void get_xy(long l, int &u, int &v) {
	u = (l >> 14);
	v = (l - (u << 14));
}

long coord_of(int u, int v) {
	long ret;
	ret = u;
	ret = ret << 14;
	ret += v;
	return ret;
}

long fill(int no_test, int nr, int nrl, int x1, int y1, long max) {	// flood-fill of line regions
	uchar thecolor;
	uchar black;
	uchar oldcolor;
	int u;
	int v;
	int ind;
	int min_x = x_global - 1;
	int min_y = y_global - 1;
	int max_x = 0;
	int max_y = 0;
	long i;
	long pos;
	long count = max;
	long count2 = 0;
	long closeby = 0;
	long *position;
	Stack<long> thestack;
	Stack<long> thestack2;
	Stack<long> thestack3;
	if (no_test) {
		count += headersize;
	}
	thecolor = i_uchar(0);
	black = i_uchar(255);
	oldcolor = theimage[nr].at<uchar>(y1, x1);
	thestack.push(coord_of(x1, y1));
	if (no_test) {
		structure[max + 8] = (long)uc_int(oldcolor);
		structure[max + 6] = numpart;
	}
	while (thestack.pop(pos)) {
		if (no_test) {
			thestack3.push(pos);
		}
		get_xy(pos, u, v);
		if (u < min_x) min_x = u;
		if (u > max_x) max_x = u;
		if (v < min_y) min_y = v;
		if (v > max_y) max_y = v;
		theimage[nr].at<uchar>(v, u) = thecolor;
		if (no_test) {
			structure[count] = pos;
		}
		count += 1;
		count2 += 1;
		if (u < x_global - 1) {
			if (theimage[nr].at<uchar>(v, u + 1) != oldcolor) {
				thestack2.push(coord_of(u + 1, v));
				closeby += 1;
			}
			if (v < y_global - 1) {
				if (theimage[nr].at<uchar>(v + 1, u + 1) != oldcolor) {
					thestack2.push(coord_of(u + 1, v + 1));
					closeby += 1;
				}
			}
		}
		if (v < y_global - 1) {
			if (theimage[nr].at<uchar>(v + 1, u) != oldcolor) {
				thestack2.push(coord_of(u, v + 1));
				closeby += 1;
			}
		}
		if (v > 0) {
			if (theimage[nr].at<uchar>(v - 1, u) == oldcolor) {
				thestack.push(coord_of(u, v - 1));
				theimage[nr].at<uchar>(v - 1, u) = black;
			}
		}
		if (v < y_global - 1) {
			if (theimage[nr].at<uchar>(v + 1, u) == oldcolor) {
				thestack.push(coord_of(u, v + 1));
				theimage[nr].at<uchar>(v + 1, u) = black;
			}
		}
		if (u > 0) {
			if (theimage[nr].at<uchar>(v, u - 1) == oldcolor) {
				thestack.push(coord_of(u - 1, v));
				theimage[nr].at<uchar>(v, u - 1) = black;
			}
		}
		if (u < x_global - 1) {
			if (theimage[nr].at<uchar>(v, u + 1) == oldcolor) {
				thestack.push(coord_of(u + 1, v));
				theimage[nr].at<uchar>(v, u + 1) = black;
			}
		}
	}
	if (no_test) {
		structure[max + 2] = count2;
		structure[max] = 1;
		structure[max + 4] = coord_of(min_x, min_y);
		structure[max + 5] = coord_of(max_x, max_y);
		while (thestack3.pop(pos)) {
			get_xy(pos, u, v);
			l_layer[nrl].at<int32_t>(v, u) = (max_x - min_x)*(max_x - min_x) + (max_y - min_y)*(max_y - min_y);
		}
	}
	if (closeby) {
		position = (long*)calloc(closeby, sizeof(long));
		closeby = 0;
		while (thestack2.pop(pos)) {
			i = 0;
			ind = 1;
			while ((ind) && (i < closeby)) {
				if (pos == position[i]) {
					ind = 0;
				}
				i += 1;
			}
			if (ind) {
				position[closeby] = pos;
				if (no_test) {
					structure[count + closeby] = pos;
				}
				closeby += 1;
			}
		}
		free(position);
	}
	if (no_test) {
		structure[max + 3] = closeby;
	}
	return count + closeby;
}

void dosort(int left, int right) {	// quick sort
	int li = left;
	int re = right;
	int test = grayvalues[(li + re) / 2];
	int mhelp;
	while (li <= re) {
		while (grayvalues[li] < test) li += 1;
		while (grayvalues[re] > test) re -= 1;
		if (li <= re) {
			mhelp = grayvalues[li];
			grayvalues[li] = grayvalues[re];
			grayvalues[re] = mhelp;
			li += 1;
			re -= 1;
		}
	}
	if (left < re) dosort(left, re);
	if (li < right) dosort(li, right);
}

float myabsvalag(float d) {
	float ret = d;
	if (d < 0) ret = -d;
	return ret;
}

int householder(long nspme, long nze, float *A, float *x) {	// least-squares using Householder method
	int ret = 0;
	long ij;
	long jj;
	long ik;
	long nsp;
	long i;
	long j;
	long k;
	float Aij;
	float Ajj;
	float beta;
	float s;
	float sigma;
	float sum;
	nsp = nspme + 1;
	for (j = 0; j < nspme; j++) {
		sigma = 0.0;
		for (i = j; i < nze; i++) {
			ij = i * nsp + j;
			Aij = A[ij];
			sigma += Aij * Aij;
		}
		if (myabsvalag(sigma) < 1e-20) return 1;
		jj = j * nsp + j;
		Ajj = A[jj];
		(Ajj < 0.0) ? (s = sqrt(sigma)) : (s = -sqrt(sigma));
		x[j] = s;
		beta = ((float) 1.0) / (s*Ajj - sigma);
		A[jj] = Ajj - s;
		for (k = j + 1; k < nsp; k++) {
			sum = 0.0;
			for (i = j; i < nze; i++) {
				ij = i * nsp + j;
				ik = i * nsp + k;
				sum += A[ij] * A[ik];
			}
			sum = beta * sum;
			for (i = j; i < nze; i++) {
				ij = i * nsp + j;
				ik = i * nsp + k;
				A[ik] += A[ij] * sum;
			}
		}
	}
	for (i = nspme - 1; i >= 0; i--) {
		sum = 0;
		for (j = i + 1; j < nspme; j++) {
			ij = i * nsp + j;
			sum += A[ij] * x[j];
		}
		k = i * nsp + nspme;
		sigma = x[i];
		if (myabsvalag(sigma) < 1e-20) {
			ret = 1;
			i = -1;
		}
		else {
			x[i] = (A[k] - sum) / sigma;
		}
	}
	return ret;
}

void determine_values(double ad1, double ad2, int &u, int &v) {
	u = (int)(ad1 + 0.5);
	v = (int)(ad2 + 0.5);
}

void lined(int nr, int x1_o, int y1_o, int x2_o, int y2_o, uchar thecolor) { // line drawing
	Point pt1;
	Point pt2;
	int i;
	pt1.x = x1_o;
	pt1.y = y1_o;
	pt2.x = x2_o;
	pt2.y = y2_o;
	LineIterator it(theimage[nr], pt1, pt2, 8, false);
	for (i = 0; i < it.count; i++, ++it) {
		theimage[nr].at<uchar>(it.pos()) = thecolor;
	}	
}

int gradients(int idx_base, int idx_grad, int idx1, int idx2, int idx1ii, int idx2ii) { // compute gradients and their direction
	int i;
	int j;
	int gv;
	int gh;
	int ret = 0;
	double w1;
	double w2;
	double absval;
	double twopi = 2 * mypi;
	double wstep;
	wstep = twopi / ((double)numpart);
	for (i = 0; i < (x_global - 1); i++) {
		for (j = 0; j < (y_global - 1); j++) {
			gv = -uc_int(theimage[idx_base].at<uchar>(j, i)) - uc_int(theimage[idx_base].at<uchar>(j + 1, i)) + uc_int(theimage[idx_base].at<uchar>(j, i + 1)) + uc_int(theimage[idx_base].at<uchar>(j + 1, i + 1));
			gh = -uc_int(theimage[idx_base].at<uchar>(j, i)) + uc_int(theimage[idx_base].at<uchar>(j + 1, i)) - uc_int(theimage[idx_base].at<uchar>(j, i + 1)) + uc_int(theimage[idx_base].at<uchar>(j + 1, i + 1));
			absval = 0.5*sqrt(gv*gv + gh * gh);
			theimage[idx_grad].at<uchar>(j, i) = i_uchar((int)absval);
			if (absval < mingrad) {
				theimage[idx1].at<uchar>(j, i) = i_uchar(0);
				theimage[idx2].at<uchar>(j, i) = i_uchar(0);
				theimage[idx1ii].at<uchar>(j, i) = i_uchar(0);
				theimage[idx2ii].at<uchar>(j, i) = i_uchar(0);
			}
			else {
				w1 = atan2((double)gv, (double)gh);
				while (w1 >= mypi) {
					w1 = w1 - twopi;
				}
				while (w1 < -mypi) {
					w1 = w1 + twopi;
				}
				w2 = w1 + mypi / ((double)numpart);
				while (w2 >= mypi) {
					w2 = w2 - twopi;
				}
				while (w2 < -mypi) {
					w2 = w2 + twopi;
				}
				theimage[idx1].at<uchar>(j, i) = i_uchar(2 * partitionvalue(w1, wstep) - 1);
				theimage[idx2].at<uchar>(j, i) = i_uchar(2 * partitionvalue(w2, wstep));
				theimage[idx1ii].at<uchar>(j, i) = theimage[idx1].at<uchar>(j, i);
				theimage[idx2ii].at<uchar>(j, i) = theimage[idx2].at<uchar>(j, i);
			}
		}
	}
	return ret;
}

int line_detection_burns(int idx_base, int idx_line, int idx1, int idx2, int idx1ii, int idx2ii, int idxl1, int idxl2) { // detection of straight lines, see http://dx.doi.org/10.1109/TPAMI.1986.4767808 
	int ret = 0;
	int i;
	int j;
	int u;
	int v;
	int m;
	int wu;
	int wv;
	int gxmin;
	int gxmax;
	int gymin;
	int gymax;
	int count2;
	int indi[4];
	int contrast;
	long k;
	long l;
	long max;
	long max_g;
	long vote;
	long count;
	float gw;
	float *A2;
	float *x;
	double ad[4][2];
	double absval;
	double gxmind;
	double gxmaxd;
	double gymind;
	double gymaxd;
	double gw_med;

	max = 0;
	count = 0;
	for (i = 0; i < x_global; i++) {	// test number of line regions and max size for memory requirements
		for (j = 0; j < y_global; j++) {
			if (uc_int(theimage[idx1ii].at<uchar>(j, i)) > 0) {
				max_g = fill(0, idx1ii, idxl1, i, j, max);
				if ((max_g - max) >= minsize) {
					count += 1;
					max = max_g;
				}
			}
			if (uc_int(theimage[idx2ii].at<uchar>(j, i)) > 0) {
				max_g = fill(0, idx2ii, idxl2, i, j, max);
				if ((max_g - max) >= minsize) {
					count += 1;
					max = max_g;
				}
			}
		}
	}
	structure_set = count * headersize + max;	// allocate data structure for lines
	structure = (long*)calloc(structure_set + headersize + minsize + 1, sizeof(long));
	if (!(structure == NULL)) {
		k = 0;
		count = 0;
		for (i = 0; i < x_global; i++) {	// actually detect and store line regions
			for (j = 0; j < y_global; j++) {
				if (uc_int(theimage[idx1].at<uchar>(j, i)) > 0) {
					max_g = fill(1, idx1, idxl1, i, j, k);
					if ((max_g - k - headersize) >= minsize) {
						structure[k + 1] = 0;
						k = max_g;
						count += 1;
					}
				}
				if (uc_int(theimage[idx2].at<uchar>(j, i)) > 0) {
					max_g = fill(1, idx2, idxl2, i, j, k);
					if ((max_g - k - headersize) >= minsize) {
						structure[k + 1] = 1;
						k = max_g;
						count += 1;
					}
				}
			}
		}
		for (i = 0; i < x_global; i++) {
			for (j = 0; j < y_global; j++) {
				if ((l_layer[0].at<int32_t>(j, i)) || (l_layer[1].at<int32_t>(j, i))) {
					if (l_layer[0].at<int32_t>(j, i) < l_layer[1].at<int32_t>(j, i)) {
						theimage[idx1].at<uchar>(j, i) = i_uchar(2);
					}
					else {
						theimage[idx1].at<uchar>(j, i) = i_uchar(1);
					}
				}
			}
		}
		k = 0;
		max_g = 0;
		while (k < structure_set) {
			if (structure[k]) {
				vote = 0;
				for (l = 0; l < structure[k + 2]; l++) {
					get_xy(structure[k + headersize + l], u, v);
					if (uc_int(theimage[idx1].at<uchar>(v, u)) == ((int)(structure[k + 1] + 1))) vote += 1;
				}
				structure[k + 9] = vote;
				if (((1.0*vote) / (1.0*structure[k + 2])) < minsup) {
					structure[k] = 0;
				}
			}
			if ((structure[k + 2] + structure[k + 3]) > max_g) max_g = structure[k + 2] + structure[k + 3];
			k += headersize + structure[k + 2] + structure[k + 3];
		}
		grayvalues = (int*)calloc(max_g, sizeof(int));
		k = 0;
		while (k < structure_set) {
			if (structure[k]) {
				for (l = 0; l < structure[k + 2] + structure[k + 3]; l++) {
					get_xy(structure[k + headersize + l], u, v);
					grayvalues[l] = uc_int(theimage[idx_base].at<uchar>(v, u));
				}
				l = structure[k + 2] + structure[k + 3];
				dosort(0, l - 1);
				if (!(l >> 3)) {
					contrast = grayvalues[l - 1] - grayvalues[0];
				}
				else {
					contrast = 0;
					for (j = 0; j < (l >> 3); j++) {
						contrast = contrast + grayvalues[l - j - 1] - grayvalues[j];
					}
					contrast = contrast / (l >> 3);
				}
				structure[k + 7] = (long)contrast;
				if (contrast < mincon) {
					structure[k] = 0;
				}
			}
			k += headersize + structure[k + 2] + structure[k + 3];
		}
		free(grayvalues);
		k = 0;
		A2 = (float*)calloc((max_g << 2), sizeof(float));
		x = (float*)calloc(3, sizeof(float));
		while (k < structure_set) {
			if (structure[k]) {
				gxmin = x_global - 1;
				gxmax = 0;
				gymin = y_global - 1;
				gymax = 0;
				m = structure[k + 2] + structure[k + 3];
				gw_med = 0;
				for (l = 0; l < m; l++) {
					get_xy(structure[k + headersize + l], u, v);
					if (u < gxmin) gxmin = u;
					if (u > gxmax) gxmax = u;
					if (v < gymin) gymin = v;
					if (v > gymax) gymax = v;
					A2[(l << 2)] = (float)u;
					A2[(l << 2) + 1] = (float)v;
					A2[(l << 2) + 2] = 1.0;
					gw = (float)uc_int(theimage[idx_base].at<uchar>(v, u));
					A2[(l << 2) + 3] = gw;
					gw_med += gw;
				}
				gw_med /= ((double)m);
				gxmind = ((double)gxmin) - 0.475;
				gymind = ((double)gymin) - 0.475;
				gxmaxd = ((double)gxmax) + 0.475;
				gymaxd = ((double)gymax) + 0.475;
				householder(3, m, A2, x);
				for (i = 0; i < 3; i++) {
					structure[k + 10 + i] = (long)(10000.0*x[i]);
				}
				ad[0][0] = (gw_med - x[2] - x[1] * gymind) / x[0];
				ad[0][1] = gymind;
				ad[1][0] = (gw_med - x[2] - x[1] * gymaxd) / x[0];
				ad[1][1] = gymaxd;
				ad[2][0] = gxmind;
				ad[2][1] = (gw_med - x[2] - x[0] * gxmind) / x[1];
				ad[3][0] = gxmaxd;
				ad[3][1] = (gw_med - x[2] - x[0] * gxmaxd) / x[1];
				count2 = 0;
				for (i = 0; i < 4; i++) indi[i] = 1;
				if ((ad[0][0] < gxmind) || (ad[0][0] > gxmaxd)) indi[0] = 0;
				if ((ad[1][0] < gxmind) || (ad[1][0] > gxmaxd)) indi[1] = 0;
				if ((ad[2][1] < gymind) || (ad[2][1] > gymaxd)) indi[2] = 0;
				if ((ad[3][1] < gymind) || (ad[3][1] > gymaxd)) indi[3] = 0;
				for (i = 0; i < 4; i++) {
					if (indi[i]) count2 += 1;
				}
				if (count2 > 1) {
					i = 0;
					while (!indi[i]) i++;
					determine_values(ad[i][0], ad[i][1], u, v);
					i += 1;
					while (((!indi[i]) || ((u == ad[i][0]) && (v == ad[i][1]))) && (i < 4)) i++;
					if (i == 4) count2 = 0;
					else {
						determine_values(ad[i][0], ad[i][1], wu, wv);
						structure[k + 13] = coord_of(u, v);
						structure[k + 14] = coord_of(wu, wv);
						absval = sqrt((u - wu)*(u - wu) + (v - wv)*(v - wv));
						if (absval < minlen) count2 = 0;
					}
				}
				if (count2 < 2) {
					structure[k] = 0;
				}
			}
			k += headersize + structure[k + 2] + structure[k + 3];
		}
		free(A2);
		free(x);
		k = 0;
		while (k < structure_set) {	// draw lines in image
			if (structure[k]) {
				get_xy(structure[k + 13], u, v);
				get_xy(structure[k + 14], i, j);
				lined(idx_line, u, v, i, j, i_uchar((int)(structure[k + 8])));
			}
			k += headersize + structure[k + 2] + structure[k + 3];
		}
		result_set = 1;
	}
	return ret;
}

int write_text_result(string mfilename) {
	int u;
	int v;
	int ret = 0;
	long k;
	long nr;
	fstream textoutp;
	if (structure_set) {
		nr = 0;
		textoutp.open(mfilename, ios::out);
		textoutp << "counter\nstart point (x/y)\nend point (x/y)\ncontrast\npartition (of " << 2 * numpart << ")" << "\n#pixels in region\n\n";
		k = 0;
		while (k < structure_set) {
			if (structure[k]) {
				textoutp << nr << endl;
				nr += 1;
				get_xy(structure[k + 13], u, v);
				textoutp << u << "\t" << v << endl;
				get_xy(structure[k + 14], u, v);
				textoutp << u << "\t" << v << endl;
				textoutp << structure[k + 7] << endl;
				textoutp << structure[k + 8] << endl;
				textoutp << structure[k + 2] + structure[k + 3] << endl;
				textoutp << endl;
			}
			k += headersize + structure[k + 2] + structure[k + 3];
		}
		textoutp.close();
		textoutp.clear();
	}
	else {
		ret = 1;
	}
	return ret;
}

int main() {
	int i;
	int j;
	double angl;
	const string inname = "./Ullsteinhaus.jpg";
	const string outname = "./Ullsteinhaus_Linien.png";
	const string outtxtname = "./Ullsteinhaus_Linien.txt";
	Mat outimg;

	read_theimage(0, x_global, y_global, inname);
	gauss(3, 4);
	gradients(4, 5, 7, 8, 9, 10);
	line_detection_burns(4, 6, 7, 8, 9, 10, 0, 1);
	write_text_result(outtxtname);	// write list of detected lines as text file

	outimg = Mat::zeros(y_global, x_global, CV_8UC3);	// write result as image using OpenCV
	for (j = 0; j < y_global; j++) {
		for (i = 0; i < x_global; i++) {
			if (theimage[6].at<uchar>(j, i)) {
				angl = mypi * ((double)(theimage[6].at<uchar>(j, i))) / ((double)numpart);
				outimg.at<Vec3b>(j, i)[0] = (uchar)(1.25 + 127.25*(1 + cos(2.0*mypi / 3.0 + angl)));
				outimg.at<Vec3b>(j, i)[1] = (uchar)(1.25 + 127.25*(1 + cos(4.0*mypi / 3.0 + angl)));
				outimg.at<Vec3b>(j, i)[2] = (uchar)(1.25 + 127.25*(1 + cos(6.0*mypi / 3.0 + angl)));
			}
		}
	}
	imwrite(outname, outimg);

	if (memory_set) {	// free memory
		for (i = 0; i < layer_c_intern; i++) {
			theimage[i].release();
		}
		for (i = 0; i < layer_l_intern; i++) {
			l_layer[i].release();
		}
		memory_set = 0;
		if (structure_set) {
			free(structure);
			structure_set = 0;
		}
		result_set = 0;
		x_global = 0;
		y_global = 0;
	}
	outimg.release();
	
	return 0;
}