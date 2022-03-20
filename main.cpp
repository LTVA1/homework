#include <QCoreApplication>
#include <QDebug>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <process.h>
#include <chrono>
#include <string>
#include <time.h>
#include <pthread.h>
#include <complex>

#define NUMBER_OF_ITERATIONS 30 /*30*/
#define MAX_STEPS 30
#define NUM_OF_THREADS 2

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
	{
		L[i] = arr[l + i];
	}

	for (j = 0; j < n2; j++)
	{
		R[j] = arr[m + 1 + j];
	}

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray
	j = 0; // Initial index of second subarray
	k = l; // Initial index of merged subarray
	while (i < n1 && j < n2)
	{
		if (L[i] <= R[j])
		{
			arr[k] = L[i];
			i++;
		}

		else
		{
			arr[k] = R[j];
			j++;
		}

		k++;
	}

	/* Copy the remaining elements of L[], if there
	are any */
	while (i < n1)
	{
		arr[k] = L[i];
		i++;
		k++;
	}

	/* Copy the remaining elements of R[], if there
	are any */
	while (j < n2)
	{
		arr[k] = R[j];
		j++;
		k++;
	}

	free(L);
	free(R);
}

void merge_sort(int* a, int r, int l)
{
	if (l < r)
	{
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

void bubblesort_swap_check(int* a, int len, int unused)
{
	len++;

	//sorted is initially false
	int comps = len - 1;
	int sorted = 0;

	while (!sorted)	//comps reduces on each pass
	{
		sorted = 1; //set true for each pass

		for (int i = 0; i < comps; i++)
		{
			if (a[i] > a[i + 1])
			{
				int temp = a[i];
				a[i] = a[i + 1];
				a[i + 1] = temp;
				sorted = 0;	//not yet sorted
			}
		}	//end of each pass

		comps--;
	}
}

struct mergeParams
{
	int begin;
	int mid;
	int end;
	int* a;
	int* sorted_a;
};

struct insertionSortParams
{
	int* a;
	int start;
	int end;
};

//int unsorted[5000000];
//int sorted[5000000];

void* merge(void* args)
{
	//unpack parameters
	struct mergeParams* params = (struct mergeParams*)args;
	int begin = params->begin,
		mid = params->mid,
		end = params->end;

	int* unsorted = params->a;
	int* sorted = params->sorted_a;

	int i = begin, j = mid, tpos = begin;

	while (i < mid && j <= end)
	{
		if (unsorted[i] < unsorted[j])
		{
			sorted[tpos++] = unsorted[i++];
		}

		else
		{
			sorted[tpos++] = unsorted[j++];
		}
	}

	//still elements left over in first list. copy over
	while (i < mid)
	{
		sorted[tpos++] = unsorted[i++];
	}

	//still elements left over in first list. copy over
	while (j <= end)
	{
		sorted[tpos++] = unsorted[j++];
	}

	return NULL;
}

void* insertionSort(void* args)
{
	struct insertionSortParams* params = (struct insertionSortParams*)args;
	int start = params->start,
		end = params->end;

	int* unsorted = params->a;

	int i = start, j, itemToInsert;

	//everything to the left of i is sorted
	while (i <= end)
	{
		itemToInsert = unsorted[i]; //a must, or else unsorted[i] gets overwritten when shifting elements

		//since everything in this sequence is sorted,
		//starting from i, and going in reverse order, shift the elements to the right
		//until an element less than unsorted[i] is found
		j = i - 1;

		while (j >= start && itemToInsert < unsorted[j])
		{
			unsorted[j + 1] = unsorted[j];
			j--;
		}
		//insert the element into s[j]
		unsorted[j + 1] = itemToInsert;
		i++;
	}

	return NULL;
}

void insertion_sort_multithread(int* a, int len, int unused)
{
	len++;

	//for(int i = 0; i < len; ++i)
	//{
		//a[i] = unsorted[i];
	//}

	//define the indices of the two sublists
	int start1 = 0,
		end1 = len / 2,
		start2 = end1 + 1,
		end2 = len - 1;

	//2 sorting threads, and 1 merge thread
	pthread_t threads[3];

	int* sorted_a = (int*)malloc(sizeof(sorted_a) * len);

	//prepare sorting params and fire off sorting threads
	struct insertionSortParams sArgs[2];
	sArgs[0].start = start1;
	sArgs[0].end = end1;
	sArgs[0].a = a;
	pthread_create(&threads[0], NULL, insertionSort, &sArgs[0]); //deal with first sublist

	sArgs[1].start = start2;
	sArgs[1].end = end2;
	sArgs[1].a = a;
	pthread_create(&threads[1], NULL, insertionSort, &sArgs[1]); //deal with second sublist

	//wait for sorting threads to terminate
	pthread_join(threads[0], NULL);
	pthread_join(threads[1], NULL);

	//prepare params and fire off merging thread
	struct mergeParams mArgs;
	mArgs.begin = start1;
	mArgs.mid = start2;
	mArgs.end = end2;
	mArgs.a = a;
	mArgs.sorted_a = sorted_a;
	pthread_create(&threads[2], NULL, merge, &mArgs); //merge the sublists into sorted[]!

	//wait for merging thread to terminate
	pthread_join(threads[2], NULL);

	//copy array
	for (int i = 0; i < len; ++i)
	{
		a[i] = sorted_a[i];
	}

	free(sorted_a);
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
	"merge sort",
	"flash sort",
	"bubblesort with swap check",
	"multithreaded inserion sort",
};

/*struct ThreadData
{
	int* a;
	int len;
	ofstream* file_writer;
	int scale_factor;
	void (*sort)(int*, int, int);
	int sort_index;
};

bool threadFinished[8];
ThreadData data[8];*/

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

void plot(FILE* gnuplot_fd, const char* filename, const char* title, int window_number, bool plot_in_new_window)
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
void plot_mandel(FILE* gnuplot_fd, const char* filename, const char* title, int window_number, bool plot_in_new_window)
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
	fprintf(gnuplot_fd, "\' using 1:2 with points pointtype 5\n");

	fflush(gnuplot_fd);
}
#endif

#define max_iteration 30
#define max_row 400
#define max_column 800

int main(int argc, char *argv[])
{
	QCoreApplication a(argc, argv);

	/*for (int i = 0; i < 8; ++i)
	{
		threadFinished[i] = false;
	}*/

	ofstream bubsort("bubsort.txt", ios::out);
	ofstream qsort("qsort.txt", ios::out);
	ofstream inssort("inssort.txt", ios::out);
	ofstream selsort("selsort.txt", ios::out);
	ofstream mergesort("mergesort.txt", ios::out);

	ofstream flashsort("flashsort.txt", ios::out);
	ofstream bubsort_swap_check("bubsort_swap_check.txt", ios::out);

	ofstream inssort_mt("inssort_mt.txt", ios::out);

#ifdef PLOT_MANDELBROT
	ofstream mandel("mandelbrot.txt", ios::out);
#endif

	for (int i = 1; i < MAX_STEPS; i++)
	{
		sort_wrap(i, bubsort, 70, bubblesort, 0);
		sort_wrap(i, qsort, 2800, quicksort, 1);
		sort_wrap(i, inssort, 140, insertion_sort, 2);
		sort_wrap(i, selsort, 100, selection_sort, 3);
		sort_wrap(i, mergesort, 1200, merge_sort, 4);
		sort_wrap(i, flashsort, 5600, flash_sort, 5);
		sort_wrap(i, bubsort_swap_check, 70, bubblesort_swap_check, 6);

		sort_wrap(i, inssort_mt, 150 * NUM_OF_THREADS, insertion_sort_multithread, 7);

		cout << "Sorting arrays, " << i << " of " << MAX_STEPS <<"...\n";
	}

	cout << "\n\nWe had 2 standard libraries, 11 strings in for cycle, 7 sorting algorithms, a memory half-full of arrays and a whole galaxy of data, arrays, random integers, graphs... Not that we strictly needed all of that for the module homework, but once you started desperate attempts of obtaining advanced standing for the course, the tendency is to push it as far as you can. The only thing that really worried me was the multithreaded sort. There is nothing in the world more annoying and complex and shooting itself in the leg than several threads writing simultaneously into one data structure, and I knew we'd get into that rotten stuff pretty soon.\n\n";

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
	bubsort_swap_check.close();

	inssort_mt.close();

#ifdef PLOT_MANDELBROT
	mandel.close();
#endif

	FILE* gnuplot_fd;

	if ((gnuplot_fd = _popen("gnuplot\\bin\\gnuplot", "w")) == NULL)
	{
		fprintf(stderr, "Error opening pipe to gnuplot.\n");
		exit(1);
	}

#ifdef PLOT_MANDELBROT
	plot_mandel(gnuplot_fd, "mandelbrot.txt", "FRACTALS ARE BEAUTIFUL", 0, true);
#endif

	plot(gnuplot_fd, "bubsort.txt", "BUBBLE SORT", 1, true);
	plot(gnuplot_fd, "qsort.txt", "QUICK SORT", 2, true);
	plot(gnuplot_fd, "inssort.txt", "INSERTION SORT", 3, true);
	plot(gnuplot_fd, "selsort.txt", "SELECTION SORT", 4, true);
	plot(gnuplot_fd, "mergesort.txt", "MERGE SORT", 5, true);

	plot(gnuplot_fd, "flashsort.txt", "FLASH SORT", 6, true);
	plot(gnuplot_fd, "bubsort_swap_check.txt", "BUBBLE SORT WITH SWAP CHECK", 7, true);

	plot(gnuplot_fd, "inssort_mt.txt", "MULTITHREADED INSERTION SORT", 8, true);

	fprintf(gnuplot_fd, "set terminal windows 9\n");
	fprintf(gnuplot_fd, "set xlabel \"Number of elements\"\n set ylabel \"Time (us)\"\n");
	fprintf(gnuplot_fd, "set title \"SORTING ALGORITHMS\"\n");

	fprintf(gnuplot_fd, "plot \"bubsort.txt\" using 1:2 title \"BUBBLE SORT\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"qsort.txt\" using 1:2 title \"QUICK SORT\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"inssort.txt\" using 1:2 title \"INSERTION SORT\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"selsort.txt\" using 1:2 title \"SELECTION SORT\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"mergesort.txt\" using 1:2 title \"MERGE SORT\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"flashsort.txt\" using 1:2 title \"FLASH SORT\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"bubsort_swap_check.txt\" using 1:2 title \"BUBBLE SORT WITH SWAP CHECK\" with linespoints\n");
	fprintf(gnuplot_fd, "replot \"inssort_mt.txt\" using 1:2 title \"MULTITHREADED INSERTION SORT\" with linespoints\n");

	fflush(gnuplot_fd);

	system("pause");

	fprintf(gnuplot_fd, "exit\n");

	_pclose(gnuplot_fd);

	cout << "\nHomework_1.cpp has left the building.\n";

	return a.exec();
}