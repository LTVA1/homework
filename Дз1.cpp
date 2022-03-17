#include <iostream>
#include <Windows.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <chrono>
#include <string> 
#include <time.h>

#include <complex>

#define NUMBER_OF_ITERATIONS 25

//#define PLOT_MANDELBROT

typedef uint64_t Uint64;
typedef int64_t Sint64;
typedef uint32_t Uint32;
typedef int32_t Sint32;
typedef uint16_t Uint16;
typedef int16_t Sint16;
typedef uint8_t Uint8;
typedef int8_t Sint8;

using namespace std;

void bubblesort(int* a, int len, int unused)
{
	len++;

	for (int i = 0; i < len; ++i)
	{
		for (int j = i + 1; j < len; ++j)
		{
			if (a[i] > a[j])
			{
				int temp = a[i];
				a[i] = a[j];
				a[j] = temp;
			}
		}
	}
}

void quicksort(int* a, int last, int first)
{
	if (first < last)
	{
		int left = first, right = last, middle = a[(left + right) / 2];

		do //partition relative to pivot (middle)
		{
			while (a[left] < middle) left++;
			while (a[right] > middle) right--;

			if (left <= right)
			{
				int tmp = a[left];
				a[left] = a[right];
				a[right] = tmp;
				left++;
				right--;
			}
		} while (left <= right);

		quicksort(a, right, first);
		quicksort(a, last, left);
	}
}

void insertion_sort(int* a, int len, int unused)
{
	len++;

	for (int i = 1; i < len; i++)
	{
		int cur = a[i];
		int j = i;

		while ((j > 0) && (cur < a[j - 1]))
		{
			a[j] = a[j - 1];
			j--;
		}

		a[j] = cur;
	}
}

void selection_sort(int* a, int len, int unused)
{
	for (int i = 0; i < len; i++) 
	{
		int minIndex = i;

		for (int j = i + 1; j < len; j++) 
		{
			if (a[j] < a[minIndex]) 
			{
				minIndex = j;
			}
		}

		int temp = a[i];
		a[i] = a[minIndex];
		a[minIndex] = temp;
	}
}




// Merges two subarrays of arr[].
// First subarray is arr[l..m]
// Second subarray is arr[m+1..r]
void merge(int arr[], int l, int m, int r)
{
	int i, j, k;
	int n1 = m - l + 1;
	int n2 = r - m;

	/* create temp arrays */
	int* L = (int*)malloc(sizeof(L) * n1);
	int* R = (int*)malloc(sizeof(R) * n2);

	/* Copy data to temp arrays L[] and R[] */
	for (i = 0; i < n1; i++)
		L[i] = arr[l + i];
	for (j = 0; j < n2; j++)
		R[j] = arr[m + 1 + j];

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray
	j = 0; // Initial index of second subarray
	k = l; // Initial index of merged subarray
	while (i < n1 && j < n2) {
		if (L[i] <= R[j]) {
			arr[k] = L[i];
			i++;
		}
		else {
			arr[k] = R[j];
			j++;
		}
		k++;
	}

	/* Copy the remaining elements of L[], if there
	are any */
	while (i < n1) {
		arr[k] = L[i];
		i++;
		k++;
	}

	/* Copy the remaining elements of R[], if there
	are any */
	while (j < n2) {
		arr[k] = R[j];
		j++;
		k++;
	}

	free(L);
	free(R);
}

void merge_sort(int* a, int r, int l)
{
	if (l < r) {
		// Same as (l+r)/2, but avoids overflow for
		// large l and h
		int m = l + (r - l) / 2;

		// Sort first and second halves
		merge_sort(a, m, l);
		merge_sort(a, r, m + 1);

		merge(a, l, m, r);
	}
}

void flash_sort(int* a, int len, int unused) //1st arbitrary sort
{
	int* __L = (int*)malloc(sizeof(__L) * len);

	int m = len * 0.43;

	if (m <= 2)
	{
		m = 2;
	}

	for (int i = 0; i < m; ++i)
	{
		__L[i] = 0;
	}

	int Mx = a[0], Mn = a[0];

	for (int i = 1; i < len; ++i)
	{
		if (Mx < a[i])
		{
			Mx = a[i];
		}

		if (Mn > a[i])
		{
			Mn = a[i];
		}
	}

	if (Mx == Mn)
	{
		return;
	}

#define getK(x) 1ll * (m - 1) * (x - Mn) / (Mx - Mn)

	for (int i = 0; i < len; ++i)
	{
		++__L[getK(a[i])];
	}

	for (int i = 1; i < m; ++i)
	{
		__L[i] += __L[i - 1];
	}

	//step 2
	int count = 0;
	int i = 0;

	while (count < len) 
	{
		int k = getK(a[i]);

		while (i >= __L[k])
		{
			k = getK(a[++i]);
		}

		int z = a[i];

		while (i != __L[k]) 
		{
			k = getK(z);
			int y = a[__L[k] - 1];
			a[--__L[k]] = z;
			z = y;
			++count;
		}
	}

	//step 3
	for (int k = 1; k < m; ++k) 
	{
		for (int i = __L[k] - 2; i >= __L[k - 1]; --i)
		{
			if (a[i] > a[i + 1]) 
			{
				int t = a[i], 
				j = i;

				while (t > a[j + 1]) 
				{ 
					a[j] = a[j + 1]; ++j; 
				}

				a[j] = t;
			}
		}
	}

	free(__L);
}

void insertion_sort_multithread(int* a, int len, int unused)
{

}

void generate_array(int*& a, int len)
{
	a = (int*)malloc(sizeof(a) * len);

	srand(time(NULL));

	for (int i = 0; i < len; ++i)
	{
		a[i] = rand();
	}
}

int check_array(int* a, int len)
{
	for (int i = 0; i < len - 2; i += 2)
	{
		if (a[i + 1] < a[i] || a[i + 2] < a[i + 1])
		{
			return -1;
		}
	}

	return 0;
}

void delete_a(int* a)
{
	free(a);
}

const char* sorts[] = {
	"bubble sort",
	"quick sort",
	"insertion sort",
	"selection sort",
	"//",
	"flash sort",
};

struct ThreadData
{
	int* a;
	int len;
	ofstream* file_writer;
	int scale_factor;
	void (*sort)(int*, int);
	int sort_index;
};

bool threadFinished[8];
ThreadData data[8];

void sort_wrap(int len, ofstream& file_writer, int scale_factor, void (*sort)(int*, int, int), int sort_index)
{
	int* a = nullptr;

	Uint64 av_time = 0;

	file_writer << len * scale_factor << "\t";

	for (int i = 0; i < NUMBER_OF_ITERATIONS; ++i)
	{
		generate_array(a, len * scale_factor);

		auto begin = chrono::high_resolution_clock::now();
		
		sort(a, len * scale_factor - 1, 0);

		auto end = chrono::high_resolution_clock::now();

		int result = check_array(a, len * scale_factor);

		if (result == -1)
		{
			cout << "\n\nRuntime error: " << sorts[sort_index] << " failed.\n\n" << endl;
			system("pause");
			exit(1);
		}

		delete_a(a);
		av_time += chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
	}

	av_time /= NUMBER_OF_ITERATIONS * 1000;

	file_writer << av_time << "\n";

	//threadFinished[sort_index] = true;
	//_endthread();
}

void plot(FILE* gnuplot_fd, const char* filename, const char* title, int window_number)
{
	string s = to_string(window_number);

	fprintf(gnuplot_fd, "set terminal windows ");
	fprintf(gnuplot_fd, s.c_str());
	fprintf(gnuplot_fd, "\nset title \'");
	fprintf(gnuplot_fd, title);
	fprintf(gnuplot_fd, "\'\n");
	fprintf(gnuplot_fd, "set xlabel \"Number of a elements\"\nset ylabel \"Time (us)\"\n");
	fprintf(gnuplot_fd, "plot \'");
	fprintf(gnuplot_fd, filename);
	fprintf(gnuplot_fd, "\' using 1:2 with linespoints\n");

	fflush(gnuplot_fd);
}

#ifdef PLOT_MANDELBROT
void plot_mandel(FILE* gnuplot_fd, const char* filename, const char* title, int window_number)
{
	string s = to_string(window_number);

	fprintf(gnuplot_fd, "set terminal windows ");
	fprintf(gnuplot_fd, s.c_str());
	fprintf(gnuplot_fd, "\nset title \'");
	fprintf(gnuplot_fd, title);
	fprintf(gnuplot_fd, "\'\n");
	fprintf(gnuplot_fd, "set xlabel \"Re\"\nset ylabel \"Im\"\n");
	fprintf(gnuplot_fd, "plot \'");
	fprintf(gnuplot_fd, filename);
	fprintf(gnuplot_fd, "\' using 1:2 with points\n");

	fflush(gnuplot_fd);
}
#endif

#define max_iteration 30
#define max_row 300
#define max_column 300

int main()
{
	for (int i = 0; i < 8; ++i)
	{
		threadFinished[i] = false;
	}

	ofstream bubsort("bubsort.txt", ios::out);
	ofstream qsort("qsort.txt", ios::out);
	ofstream inssort("inssort.txt", ios::out);
	ofstream selsort("selsort.txt", ios::out);
	ofstream mergesort("mergesort.txt", ios::out);
	ofstream flashsort("flashsort.txt", ios::out);

#ifdef PLOT_MANDELBROT
	ofstream mandel("mandelbrot.txt", ios::out);
#endif

	for (int i = 1; i < 31; i++)
	{
		sort_wrap(i, bubsort, 70, bubblesort, 0);
		sort_wrap(i, qsort, 2600, quicksort, 1);
		sort_wrap(i, inssort, 140, insertion_sort, 2);
		sort_wrap(i, selsort, 110, selection_sort, 3);
		sort_wrap(i, mergesort, 500, merge_sort, 4);
		sort_wrap(i, flashsort, 4000, flash_sort, 5);
	}

#ifdef PLOT_MANDELBROT
	for (int row = 0; row < max_row; ++row) //a little beautiful thing
	{
		for (int column = 0; column < max_column; ++column) 
		{
			complex<double> z, c = 
			{
				(double)column * 2 / max_column - 1.5f,
				(double)row * 2 / max_row - 1
			};

			int iteration = 0;

			while (abs(z) < 2 && ++iteration < max_iteration)
			{
				z = pow(z, 2) + c;
			}

			if (iteration == max_iteration)
			{
				mandel << column << "\t" << row << "\n";
			}
		}
	}
#endif

	bubsort.close();
	qsort.close();
	inssort.close();
	selsort.close();
	mergesort.close();

	flashsort.close();

#ifdef PLOT_MANDELBROT
	mandel.close();
#endif

	FILE* gnuplot_fd;

	if ((gnuplot_fd = _popen("gnuplot\\bin\\gnuplot.exe", "w")) == NULL) //if ((gnuplot_fd = _popen("gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "Error opening pipe to gnuplot.\n");
		exit(1);
	}

#ifdef PLOT_MANDELBROT
	plot_mandel(gnuplot_fd, "mandelbrot.txt", "FRACTALS ARE BEAUTIFUL", 0);
#endif

	plot(gnuplot_fd, "bubsort.txt", "BUBBLE SORT", 1);
	plot(gnuplot_fd, "qsort.txt", "QUICK SORT", 2);
	plot(gnuplot_fd, "inssort.txt", "INSERTION SORT", 3);
	plot(gnuplot_fd, "selsort.txt", "SELECTION SORT", 4);
	plot(gnuplot_fd, "mergesort.txt", "MERGE SORT", 5);

	plot(gnuplot_fd, "flashsort.txt", "FLASH SORT", 6);

	system("pause");

	fprintf(gnuplot_fd, "exit\n");

	_pclose(gnuplot_fd);

	return 0;
}