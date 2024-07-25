#include <iostream>
#include<mpi.h>
#include<cstdlib>
using namespace std;
 
void binarySearch(int* arr, int start, int end, int key, int rank)
{
	while (start <= end)
	{
		int mid = (start + end) / 2;
		if (arr[mid] == key)
		{
			cout << "Element is  Found by processor " << rank <<  " .\n";
			return;
		}
		else if (arr[mid] < key)
		{
			start = mid + 1;
		}
		else
		{
			end = mid - 1;
		}
	}
}
 
void merge(int array[], int const left, int const mid,
	int const right)
{
	int const subArrayOne = mid - left + 1;
	int const subArrayTwo = right - mid;
 
	// Create temp arrays
	auto* leftArray = new int[subArrayOne],
		* rightArray = new int[subArrayTwo];
 
	// Copy data to temp arrays leftArray[] and rightArray[]
	for (auto i = 0; i < subArrayOne; i++)
		leftArray[i] = array[left + i];
	for (auto j = 0; j < subArrayTwo; j++)
		rightArray[j] = array[mid + 1 + j];
 
	auto indexOfSubArrayOne = 0, indexOfSubArrayTwo = 0;
	int indexOfMergedArray = left;
 
	// Merge the temp arrays back into array[left..right]
	while (indexOfSubArrayOne < subArrayOne
&& indexOfSubArrayTwo < subArrayTwo) {
		if (leftArray[indexOfSubArrayOne]
<= rightArray[indexOfSubArrayTwo]) {
			array[indexOfMergedArray]
				= leftArray[indexOfSubArrayOne];
			indexOfSubArrayOne++;
		}
		else {
			array[indexOfMergedArray]
				= rightArray[indexOfSubArrayTwo];
			indexOfSubArrayTwo++;
		}
		indexOfMergedArray++;
	}
 
	// Copy the remaining elements of
	// left[], if there are any
	while (indexOfSubArrayOne < subArrayOne) {
		array[indexOfMergedArray]
			= leftArray[indexOfSubArrayOne];
		indexOfSubArrayOne++;
		indexOfMergedArray++;
	}
 
	// Copy the remaining elements of
	// right[], if there are any
	while (indexOfSubArrayTwo < subArrayTwo) {
		array[indexOfMergedArray]
			= rightArray[indexOfSubArrayTwo];
		indexOfSubArrayTwo++;
		indexOfMergedArray++;
	}
	delete[] leftArray;
	delete[] rightArray;
}
 
void mergeSort(int array[], int const begin, int const end)
{
	if (begin >= end)
		return;
 
	int mid = begin + (end - begin) / 2;
	mergeSort(array, begin, mid);
	mergeSort(array, mid + 1, end);
	merge(array, begin, mid, end);
}
 
int main(int argc, char** argv) {
	int n = 4000;
	int* arr = new int[n];
 
	for (int i = 0; i < n; i++)
	{
		arr[i] = rand() % n;
	}
 
	int key; //element to search
	MPI_Init(&argc, &argv);
	int rank, size;double parallel = 0,parallel2,sequential;
 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
 
	if (rank == 0) {
		cout << "Enter element to be found:";
		cin >> key;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&key, 1, MPI_INT, 0, MPI_COMM_WORLD);
 
	int blocks = 6; //number of processors in paralell
	int blockSize = n / blocks;
 
 
	if (rank == 0)
	{
		double start = MPI_Wtime();
		mergeSort(arr, rank * blockSize, (rank + 1) * blockSize - 1);
		binarySearch(arr, rank * blockSize, (rank + 1) * blockSize - 1, key, rank);
		double end = MPI_Wtime();
		parallel = (end - start) * 1000;
		cout << "Execution time of Processor " << 0 << " is " << parallel << endl;
		for (int i = 1; i < 6; i++) {
			MPI_Recv(&parallel2, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			cout << "Execution time of Processor " << i << " is " << parallel2 << endl;
			parallel = max(parallel, parallel2);
		}
		MPI_Recv(&sequential, 1, MPI_DOUBLE, 6, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		cout << "Execution time of Processor " << 6 << "(sequential) is " << sequential << endl;
		cout << endl;
		cout << "Speed up :" << sequential / parallel<<endl;
		cout << "Efficiency:" << sequential / parallel / 6<<endl;
	}
	else if (rank == 1)
	{
		double start = MPI_Wtime();
		mergeSort(arr, rank * blockSize, (rank + 1) * blockSize - 1);
		binarySearch(arr, rank * blockSize, (rank + 1) * blockSize - 1, key, rank);
		double end = MPI_Wtime();
		parallel2 = (end - start) * 1000;
		MPI_Send(&parallel2,1,MPI_DOUBLE,0,rank,MPI_COMM_WORLD);
	}
	else if (rank == 2)
	{
		double start = MPI_Wtime();
		mergeSort(arr, rank * blockSize, (rank + 1) * blockSize - 1);
		binarySearch(arr, rank * blockSize, (rank + 1) * blockSize - 1, key, rank);
		double end = MPI_Wtime();
		parallel2 = (end - start) * 1000;
		MPI_Send(&parallel2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	else if (rank == 3)
	{
		double start = MPI_Wtime();
		mergeSort(arr, rank * blockSize, (rank + 1) * blockSize - 1);
		binarySearch(arr, rank * blockSize, (rank + 1) * blockSize - 1, key, rank);
		double end = MPI_Wtime();
		parallel2 = (end - start) * 1000;
		MPI_Send(&parallel2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	else if (rank == 4)
	{
		double start = MPI_Wtime();
		mergeSort(arr, rank * blockSize, (rank + 1) * blockSize - 1);
		binarySearch(arr, rank * blockSize, (rank + 1) * blockSize - 1, key, rank);
		double end = MPI_Wtime();
		parallel2 = (end - start) * 1000;
		MPI_Send(&parallel2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	else if (rank == 5)
	{
		double start = MPI_Wtime();
		mergeSort(arr, rank * blockSize, (rank + 1) * blockSize - 1);
		binarySearch(arr, rank * blockSize, (rank + 1) * blockSize - 1, key, rank);
		double end = MPI_Wtime();
		parallel2 = (end - start) * 1000;
		MPI_Send(&parallel2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	else if (rank == 6)
	{
		double start = MPI_Wtime();
		mergeSort(arr, 0,4000);
		binarySearch(arr, 0,4000, key, rank);
		double end = MPI_Wtime();
		parallel2 = (end - start) * 1000;
		MPI_Send(&parallel2, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	}
	MPI_Finalize();
 
	return 0;
}
