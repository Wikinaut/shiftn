import math
import os
import curses
import numpy as np

# implementation of line detection as described in http://dx.doi.org/10.1109/TPAMI.1986.4767808
# Marcus Hebel

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <fstream>
#include <iostream>
#include "mystack.h"

#using namespace std
#using namespace cv

layer_c_intern = 11
layer_l_intern = 2
mypi = 3.14159265358979323846264338328

# line extraction parameter

numpart = 10	# number of partitions of 360 degrees
minsize = 46	# min region size for line (approx. double line lenghth)
minlen = 23	# min line length
mincon = 3	# minimum contrast of edges
mingrad = 1.5	# minimum gradient of edges
minsup = 0.6	# minimum portion of pixels supporting the linergb = ( 30,59,11 ); // rgb portions to convert to gray-scale

#Mat theimage[layer_c_intern]	# data structure to store the image channels (r, g, b, gray, low-pass, gradient, lines, tmp1, tmp2, tmp3, tmp4)
#Mat l_layer[layer_l_intern]	# data structure for calculations

theimage = np.array()
l_layer = np.array()

# data structure to store the list of detected lines
structure_set = 0

memory_set = 0
result_set = 0
headersize = 15
x_global = 0
y_global = 0

def i_uchar(i):
	if i > 255: i = 255
	if i < 0: i = 0
	return chr(i)

def read_theimage(id, x, y, mfilename):
	global theimage
	
	# read image using OpenCV

	cvimg = np.array()

	u = 0
	v = 0
	ret = 0
	c = ""

	cvimg = imread(mfilename)
	u = cvimg.cols
	v = cvimg.rows

	if (u == 0) or (v == 0) or (u > 65535) or (v > 65535): ret = 1

	if not ret:
		for i in range(layer_c_intern):
			theimage[i] = np.zeros(v, u, dtype=int)
		for i in range(layer_l_intern):
			l_layer[i] = np.zeros(v, u, dtype=float)

	if not ret:
		memory_set = 1
		for j in range(v):
			for i in range(u):
				data = 0
				#c = cvimg.at<Vec3b>(j, i)[0] # ???????????????????????
				c = cvimg(j, i)[0]
				#theimage[id].at<uchar>(j, i) = c
				theimage[id](j, i) = c
				data += rgb[2] * ord(c)
				c = cvimg.at<Vec3b>(j, i)[1]
				theimage[id + 1].at<uchar>(j, i) = c
				data += rgb[1] * ord(c)
				c = cvimg.at<Vec3b>(j, i)[2]
				theimage[id + 2].at<uchar>(j, i) = c
				data += rgb[0] * ord(c)
				data /= (rgb[0] + rgb[1] + rgb[2])
				theimage[id + 3].at<uchar>(j, i) = i_uchar(data)
	x = u
	y = v
	cvimg.release()
	return ret


def gauss(id_old, id_new):
	# low-pass filter of grayscale image

	for j in range(y_global):
		for i in range(x_global):
			thesum = 0
			thecount = 0
			for k in range(-1, 1+1):
				for l in range(-1, 1+1):
					if (i + k >= 0) and (i + k < x_global) and (j + l >= 0) and (j + l < y_global):
						themultiplicator = 3 - k * k - l * l
						if (k == 0) and (l == 0): themultiplicator += 1
						thesum += themultiplicator * ord(theimage[id_old].at<uchar>(j + l, i + k))
						thecount += themultiplicator
			theimage[id_new].at<uchar>(j, i) = i_uchar(int((thesum / thecount)))

	ret = 0

	wdown = -mypi
	i = 1
	while i <= numpart:
		if (w >= wdown) and (w < (wdown + wsteph)):
			ret = i
			i = numpart
		wdown += wsteph
		i += 1
	return ret

def get_xy(l, u, v):
	u = (l >> 14)
	v = (l - (u << 14))

	ret = u
	ret = ret << 14
	ret += v
	return ret

def fill(no_test, nr, nrl, x1, y1, max):
	# flood-fill of line regions
	thecolor=""
	black=""
	oldcolor=""

	min_x = x_global - 1
	min_y = y_global - 1
	max_x = 0
	max_y = 0

	count = max
	count2 = 0
	closeby = 0

	Stack<long> thestack
	Stack<long> thestack2
	Stack<long> thestack3
	if no_test:
		count += headersize
	thecolor = i_uchar(0)
	black = i_uchar(255)
	oldcolor = theimage[nr].at<uchar>(y1, x1)
	thestack.push(coord_of(x1, y1))
	if no_test:
		structure[max + 8] = int(oldcolor)
		structure[max + 6] = numpart
	while thestack.pop(pos):
		if no_test:
			thestack3.push(pos)
		get_xy(pos, u, v)
		if u < min_x: min_x = u
		if u > max_x: max_x = u
		if v < min_y: min_y = v
		if v > max_y: max_y = v
		theimage[nr].at<uchar>(v, u) = thecolor
		if no_test:
			structure[count] = pos
		count += 1
		count2 += 1
		if u < x_global - 1:
			if theimage[nr].at<uchar>(v, u + 1) != oldcolor:
				thestack2.push(coord_of(u + 1, v))
				closeby += 1
			if v < y_global - 1:
				if theimage[nr].at<uchar>(v + 1, u + 1) != oldcolor:
					thestack2.push(coord_of(u + 1, v + 1))
					closeby += 1
		if v < y_global - 1:
			if theimage[nr].at<uchar>(v + 1, u) != oldcolor:
				thestack2.push(coord_of(u, v + 1))
				closeby += 1
		if v > 0:
			if theimage[nr].at<uchar>(v - 1, u) == oldcolor:
				thestack.push(coord_of(u, v - 1))
				theimage[nr].at<uchar>(v - 1, u) = black
		if v < y_global - 1:
			if theimage[nr].at<uchar>(v + 1, u) == oldcolor:
				thestack.push(coord_of(u, v + 1))
				theimage[nr].at<uchar>(v + 1, u) = black
		if u > 0:
			if theimage[nr].at<uchar>(v, u - 1) == oldcolor:
				thestack.push(coord_of(u - 1, v))
				theimage[nr].at<uchar>(v, u - 1) = black
		if u < x_global - 1:
			if theimage[nr].at<uchar>(v, u + 1) == oldcolor:
				thestack.push(coord_of(u + 1, v))
				theimage[nr].at<uchar>(v, u + 1) = black
	if no_test:
		structure[max + 2] = count2
		structure[max] = 1
		structure[max + 4] = coord_of(min_x, min_y)
		structure[max + 5] = coord_of(max_x, max_y)
		while thestack3.pop(pos):
			get_xy(pos, u, v)
			l_layer[nrl].at<int32_t>(v, u) = (max_x - min_x)*(max_x - min_x) + (max_y - min_y)*(max_y - min_y)
	if closeby:
		##### position = ()calloc(closeby, sizeof(long)); # ?????
		closeby = 0
		while thestack2.pop(pos):
			i = 0
			ind = 1
			while (ind) and (i < closeby):
				if pos == position[i]:
					ind = 0
				i += 1
			if ind:
				position[closeby] = pos
				if no_test:
					structure[count + closeby] = pos
				closeby += 1
		free(position)
	if no_test:
		structure[max + 3] = closeby
	return count + closeby

def dosort(left, right):
	# quick sort
	li = left
	re = right
	test = grayvalues[(li + re) / 2]

	while li <= re:
		while grayvalues[li] < test: li += 1
		while grayvalues[re] > test: re -= 1
		if li <= re:
			mhelp = grayvalues[li]
			grayvalues[li] = grayvalues[re]
			grayvalues[re] = mhelp
			li += 1
			re -= 1
	if left < re: dosort(left, re)
	if li < right: dosort(li, right)

def myabsvalag(d):
	ret = d
	if d < 0: ret = -d
	return ret

def householder(nspme, nze, A, x):
	# least-squares using Householder method
	ret = 0

	nsp = nspme + 1
	for j in range(nspme):
		sigma = 0.0
		for i in range(j, nze):
			ij = i * nsp + j
			Aij = A[ij]
			sigma += Aij * Aij
		if myabsvalag(sigma) < 1e-20: return 1
		jj = j * nsp + j
		Ajj = A[jj]
		s = math.sqrt(sigma) if (Ajj < 0.0) else -math.sqrt(sigma)
		x[j] = s
		beta = 1.0 / (s*Ajj - sigma)
		A[jj] = Ajj - s
		for k in range(j + 1, nsp):
			sum = 0.0
			for i in range(j, nze):
				ij = i * nsp + j
				ik = i * nsp + k
				sum += A[ij] * A[ik]
			sum = beta * sum
			for i in range(j, nze):
				ij = i * nsp + j
				ik = i * nsp + k
				A[ik] += A[ij] * sum

	# for (i = nspme - 1; i >= 0; i--):
	i = nspme-1
	while i >= 0:
		sum = 0
		for j in range(i + 1, nspme):
			ij = i * nsp + j
			sum += A[ij] * x[j]
		k = i * nsp + nspme
		sigma = x[i]
		if myabsvalag(sigma) < 1e-20:
			ret = 1
			i = -1
		else:
			x[i] = (A[k] - sum) / sigma
		i = i-1

	return ret

def determine_values(ad1, ad2, u, v):
	u = (int)(ad1 + 0.5)
	v = (int)(ad2 + 0.5)

def lined(nr, x1_o, y1_o, x2_o, y2_o, thecolor):
	# line drawing
	pt1 = {'x': x1_o, 'y': y1_o}
	pt2 = {'x': x2_o, 'y': y2_o}
	it = cv.LineIterator(theimage[nr], pt1, pt2, 8, False)
	
	# for (i = 0; i < it.count; i++, ++it): # ?????????????????
	i = 0
	while (i < it.count):
		theimage[nr].at<uchar>(it.pos()) = thecolor
		i = i + 1


def gradients(idx_base, idx_grad, idx1, idx2, idx1ii, idx2ii):
	# compute gradients and their direction

	ret = 0
	twopi = 2 * mypi

	# wstep = twopi / ((double)numpart); # ???????????????????
	wstep = twopi / numpart
	for i in range((x_global - 1)):
		for j in range((y_global - 1)):
			gv = -ord(theimage[idx_base].at<uchar>(j, i)) - ord(theimage[idx_base].at<uchar>(j + 1, i)) + ord(theimage[idx_base].at<uchar>(j, i + 1)) + ord(theimage[idx_base].at<uchar>(j + 1, i + 1))
			gh = -ord(theimage[idx_base].at<uchar>(j, i)) + ord(theimage[idx_base].at<uchar>(j + 1, i)) - ord(theimage[idx_base].at<uchar>(j, i + 1)) + ord(theimage[idx_base].at<uchar>(j + 1, i + 1))
			absval = 0.5*math.sqrt(gv*gv + gh * gh)
			theimage[idx_grad].at<uchar>(j, i) = i_uchar(absval)
			if absval < mingrad:
				theimage[idx1].at<uchar>(j, i) = i_uchar(0)
				theimage[idx2].at<uchar>(j, i) = i_uchar(0)
				theimage[idx1ii].at<uchar>(j, i) = i_uchar(0)
				theimage[idx2ii].at<uchar>(j, i) = i_uchar(0)
			else:
				w1 = math.atan2(gv, gh)
				while w1 >= mypi:
					w1 = w1 - twopi
				while w1 < -mypi:
					w1 = w1 + twopi
				w2 = w1 + mypi / numpart
				while w2 >= mypi:
					w2 = w2 - twopi
				while w2 < -mypi:
					w2 = w2 + twopi
				theimage[idx1].at<uchar>(j, i) = i_uchar(2 * partitionvalue(w1, wstep) - 1)
				theimage[idx2].at<uchar>(j, i) = i_uchar(2 * partitionvalue(w2, wstep))
				theimage[idx1ii].at<uchar>(j, i) = theimage[idx1].at<uchar>(j, i)
				theimage[idx2ii].at<uchar>(j, i) = theimage[idx2].at<uchar>(j, i)
	return ret

def line_detection_burns(idx_base, idx_line, idx1, idx2, idx1ii, idx2ii, idxl1, idxl2):
	# detection of straight lines, see http://dx.doi.org/10.1109/TPAMI.1986.4767808
	ret = 0

	max = 0
	count = 0
	for i in range(x_global):	# test number of line regions and max size for memory requirements
		for j in range(y_global):
			if ord(theimage[idx1ii].at<uchar>(j, i)) > 0:
				max_g = fill(0, idx1ii, idxl1, i, j, max)
				if (max_g - max) >= minsize:
					count += 1
					max = max_g
			if ord(theimage[idx2ii].at<uchar>(j, i)) > 0:
				max_g = fill(0, idx2ii, idxl2, i, j, max)
				if (max_g - max) >= minsize:
					count += 1
					max = max_g

	structure_set = count * headersize + max	# allocate data structure for lines
	# ???????????? structure = ()calloc(structure_set + headersize + minsize + 1, sizeof(long)); # ?????????????????????
	if not (structure == NULL):
		k = 0
		count = 0
		for i in range(x_global):	# actually detect and store line regions
			for j in range(y_global):
				if ord(theimage[idx1].at<uchar>(j, i)) > 0:
					max_g = fill(1, idx1, idxl1, i, j, k)
					if (max_g - k - headersize) >= minsize:
						structure[k + 1] = 0
						k = max_g
						count += 1
				if ord(theimage[idx2].at<uchar>(j, i)) > 0:
					max_g = fill(1, idx2, idxl2, i, j, k)
					if (max_g - k - headersize) >= minsize:
						structure[k + 1] = 1
						k = max_g
						count += 1
		for i in range(x_global):
			for j in range(y_global):
				if (l_layer[0].at<int32_t>(j, i)) or (l_layer[1].at<int32_t>(j, i)):
					if l_layer[0].at<int32_t>(j, i) < l_layer[1].at<int32_t>(j, i):
						theimage[idx1].at<uchar>(j, i) = i_uchar(2)
					else:
						theimage[idx1].at<uchar>(j, i) = i_uchar(1)
		k = 0
		max_g = 0
		while k < structure_set:
			if structure[k]:
				vote = 0
				for l in range(structure[k + 2]):
					get_xy(structure[k + headersize + l], u, v)
					if ord(theimage[idx1].at<uchar>(v, u)) == ((int)(structure[k + 1] + 1)): vote += 1
				structure[k + 9] = vote
				if ((1.0*vote) / (1.0*structure[k + 2])) < minsup:
					structure[k] = 0
			if (structure[k + 2] + structure[k + 3]) > max_g:
				max_g = structure[k + 2] + structure[k + 3]
			k += headersize + structure[k + 2] + structure[k + 3]

		# ????? grayvalues = ()calloc(max_g, sizeof(int)); # ??????????

		k = 0
		while k < structure_set:
			if structure[k]:
				for l in range(structure[k + 2] + structure[k + 3]):
					get_xy(structure[k + headersize + l], u, v)
					grayvalues[l] = ord(theimage[idx_base].at<uchar>(v, u))
				l = structure[k + 2] + structure[k + 3]
				dosort(0, l - 1)
				if not (l >> 3):
					contrast = grayvalues[l - 1] - grayvalues[0]
				else:
					contrast = 0
					for j in range((l >> 3)):
						contrast = contrast + grayvalues[l - j - 1] - grayvalues[j]
					contrast = contrast / (l >> 3)
				structure[k + 7] = contrast
				if contrast < mincon:
					structure[k] = 0
			k += headersize + structure[k + 2] + structure[k + 3]
		free(grayvalues)
		k = 0
		# ????? A2 = ()calloc((max_g << 2), sizeof(float)); # ???????????
		# ?????? x = ()calloc(3, sizeof(float)); #??????
		while k < structure_set:
			if structure[k]:
				gxmin = x_global - 1
				gxmax = 0
				gymin = y_global - 1
				gymax = 0
				m = structure[k + 2] + structure[k + 3]
				gw_med = 0
				for l in range(m):
					get_xy(structure[k + headersize + l], u, v)
					if u < gxmin: gxmin = u
					if u > gxmax: gxmax = u
					if v < gymin: gymin = v
					if v > gymax: gymax = v
					A2[(l << 2)] = float(u)
					A2[(l << 2) + 1] = float(v)
					A2[(l << 2) + 2] = 1.0
					gw = float(ord(theimage[idx_base].at<uchar>(v, u)))
					A2[(l << 2) + 3] = gw
					gw_med += gw
				gw_med /= double(m)
				gxmind = gxmin - 0.475
				gymind = gymin - 0.475
				gxmaxd = gxmax + 0.475
				gymaxd = gymax + 0.475
				householder(3, m, A2, x)
				for i in range(3):
					structure[k + 10 + i] = 10000.0*x[i]
				ad[0][0] = (gw_med - x[2] - x[1] * gymind) / x[0]
				ad[0][1] = gymind
				ad[1][0] = (gw_med - x[2] - x[1] * gymaxd) / x[0]
				ad[1][1] = gymaxd
				ad[2][0] = gxmind
				ad[2][1] = (gw_med - x[2] - x[0] * gxmind) / x[1]
				ad[3][0] = gxmaxd
				ad[3][1] = (gw_med - x[2] - x[0] * gxmaxd) / x[1]
				count2 = 0

				for i in range(0, 4): indi[i] = 1
				if (ad[0][0] < gxmind) or (ad[0][0] > gxmaxd): indi[0] = 0
				if (ad[1][0] < gxmind) or (ad[1][0] > gxmaxd): indi[1] = 0
				if (ad[2][1] < gymind) or (ad[2][1] > gymaxd): indi[2] = 0
				if (ad[3][1] < gymind) or (ad[3][1] > gymaxd): indi[3] = 0
				for i in range(4):
					if indi[i]: count2 += 1
				if count2 > 1:
					i = 0
					while not indi[i]: i += 1
					determine_values(ad[i][0], ad[i][1], u, v)
					i += 1
					while ((not indi[i]) or ((u == ad[i][0]) and (v == ad[i][1]))) and (i < 4): i += 1
					if i == 4: count2 = 0
					else:
						determine_values(ad[i][0], ad[i][1], wu, wv)
						structure[k + 13] = coord_of(u, v)
						structure[k + 14] = coord_of(wu, wv)
						absval = math.sqrt((u - wu)*(u - wu) + (v - wv)*(v - wv))
						if absval < minlen: count2 = 0
				if count2 < 2:
					structure[k] = 0
			k += headersize + structure[k + 2] + structure[k + 3]
		free(A2)
		free(x)
		k = 0
		while k < structure_set:	# draw lines in image
			if structure[k]:
				get_xy(structure[k + 13], u, v)
				get_xy(structure[k + 14], i, j)
				lined(idx_line, u, v, i, j, i_uchar((int)(structure[k + 8])))
			k += headersize + structure[k + 2] + structure[k + 3]
		result_set = 1
	return ret


def main():
	ret = 0

#	fstream textoutp
	if structure_set:
		nr = 0
#		textoutp.os.open(mfilename, ios::out)
#		textoutp << "counter\nstart point (x/y)\nend point (x/y)\ncontrast\npartition (of " << 2 * numpart << ")" << "\n#pixels in region\n\n"
		k = 0
		while k < structure_set:
			if structure[k]:
				textoutp << nr << endl
				nr += 1
				get_xy(structure[k + 13], u, v)
				textoutp << u << "\t" << v << endl
				get_xy(structure[k + 14], u, v)
				textoutp << u << "\t" << v << endl
				textoutp << structure[k + 7] << endl
				textoutp << structure[k + 8] << endl
				textoutp << structure[k + 2] + structure[k + 3] << endl
				textoutp << endl
			k += headersize + structure[k + 2] + structure[k + 3]
		textoutp.os.close()
		textoutp.stdscr.clear()
	else:
		ret = 1
	return ret

	inname = "./Ullsteinhaus.jpg"
	outname = "./Ullsteinhaus_Linien.png"
	outtxtname = "./Ullsteinhaus_Linien.txt"
	outimg = np.array()

	read_theimage(0, x_global, y_global, inname)
	gauss(3, 4)
	gradients(4, 5, 7, 8, 9, 10)
	line_detection_burns(4, 6, 7, 8, 9, 10, 0, 1)
	write_text_result(outtxtname)	# write list of detected lines as text file

	outimg = np.zeros(y_global, x_global, dtype=cv.CV_8UC3)	# write result as image using OpenCV
	for j in range(y_global):
		for i in range(x_global):
			if theimage[6].at<uchar>(j, i):
				angl = mypi * (float(theimage[6].at<uchar>(j, i))) / float(numpart)
				outimg.at<Vec3b>(j, i)[0] = (uchar)(1.25 + 127.25*(1 + math.cos(2.0*mypi / 3.0 + angl)))
				outimg.at<Vec3b>(j, i)[1] = (uchar)(1.25 + 127.25*(1 + math.cos(4.0*mypi / 3.0 + angl)))
				outimg.at<Vec3b>(j, i)[2] = (uchar)(1.25 + 127.25*(1 + math.cos(6.0*mypi / 3.0 + angl)))
	imwrite(outname, outimg)

	if memory_set:	# free memory
		for i in range(layer_c_intern):
			theimage[i].release()
		for i in range(layer_l_intern):
			l_layer[i].release()
		memory_set = 0
		if structure_set:
			free(structure)
			structure_set = 0
		result_set = 0
		x_global = 0
		y_global = 0
	outimg.release()

	return 0

if __name__ == "__main__":
    sys.exit(main())
