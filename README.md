# Start Up Code for Polynomial Class
## Purpose
The goal of this lab is to practice dynamic memory with a simple class. Additionally you will have the opportunity to practice command line compilation and memory leak checking using valgrind.

## Design
This Polynomial class will store the coefficients of a one variable polynomial  in a dynamically created array. This class has two data members:
```c++
	size_t _degree;
	float* _coefficients;
```
It is straightforward to see what they represent.
For example:
![equation](https://latex.codecogs.com/gif.latex?7x%5E4-2x%5E2&plus;0.5x&plus;2.23)
would be:
```c++
	_degree = 4
	_coefficients => {2.23, 0.5, -2.0, 0.0, 7}
```
Notice that the array size is `_degree + 1`.

## What is given?
You will be provided with the code of the constructors, the ToString, Minus, Assignment Operator, Equals, Read and Write  methods
the rest of the methods you will implement yourself.
Additionally, you will be provided with the unit tests to test your code.

## What is expected?
- Your code ***must*** compile.
- Your grade will be *almost* the percentage of your passed unit tests.
- Your code needs to show good programming practices: indentation, meaningful variable names, identifiers convention (CamelCase for functions, 
camelCase for variables, _camelCase for data members, SNAKE_CASE for constants), header comments, correct file names, etc. Failure to 
code appropriate will result in strong points penalization.
- You will need to provide implementation of the following methods:
  1. Destructor
  2. Add
  3. Subtract
  4. Multiply
  5. Derive
  6. Evaluate
- The following methods will count toward extra credit, some code will be provided to guarantee compilation:
  1. Divide
  2. Integrate

## Brief description of the methods
- `Polynomial(size_t degree)`. GIVEN. Creates a polynomial of degree `degree`, dynamically allocates an array of size `degree + 1` and sets the coefficients to 0.0
- `Polynomial(size_t degree, const float* coefficients)`. GIVEN. Creates a polynomial of degree `degree` and with the values of the coefficients in the `coefficients` array
- `Polynomial(const Polynomial& other)`. GIVEN. Copy constructor, creates a deep copy of the `other` polynomial.
- `~Polynomial`. TODO. Destructor, deallocates the allocated memory.
- `const Polynomial Sum(const Polynomial& rhs)const`. TODO. Returns a polynomial that will contain the sum of `this` with `rhs`. Notice that this operation does not change `this`.
- `const Polynomial Subtract(const Polynomial& rhs)const`. TODO. Returns a polynomial that will contain the subtraction of `this` with `rhs`. Notice that this operation does not change `this`.
- `const Polynomial Minus()const`. TODO. Returns a polynomial that will contain the additive inverse of `this`. For instance if the polynomial is p = 2x^2 - 2, then p.Minus() would return -2x^2 + 2. Notice that this operation does not change `this`.
- `const Polynomial Multiply(const Polynomial& rhs)const`. TODO. Returns a polynomial that will contain the multiplication of `this` with `rhs`. Notice that this operation does not change `this`.
- `const Polynomial Divide(const Polynomial& rhs)const`. TODO. EXTRA CREDIT. Returns a polynomial that will contain the division of `this` with `rhs`. Notice that this operation does not change `this`. The reminder is lost if there was one.
- `const Polynomial Derive()const`. TODO. Returns a polynomial that will contain the derivative of `this`. Notice that this operation does not change `this`.
- `float Evaluate(float x)const`. TODO. Returns the value of evaluating the polynomial with value `x`. Notice that this operation does not change `this`.
- `float Integrate(float a, float b)const`. TODO. EXTRA CREDIT. Returns the value of the definite integral evaluated from `a` to `b`. Notice that this operation does not change `this`.
- `const Polynomial& operator=(const Polynomial& rhs)`. GIVEN. Assigns the polynomial `rhs` to `this`.
- `string ToString()const`. GIVEN. Returns a `string` representation of the polynomial.
- `bool Equals(const Polynomial& rhs)const`. GIVEN. Returns `true` if `this` is equal to `rhs`.
- `ostream& Write(ostream& output)const`. GIVEN. Writes a polynomial from the `output` stream.
- `istream& Read(istream& input)`. GIVEN. Reads a polynomial from the `input` stream.

## How to "make", run, and memory check your program
For this assignment you are expected to compile, run and do memory check of your code in the command line. Below you will find some instructions to do so. **All the instructions assume that you are located in your project directory**.

### Compiling and running your Code
This assignment includes a `makefile`. This is a file that has the instructions on how to compile your project. The `makefile` has three options:
- Default `make`, this one would compile `polynomial.cpp` and `main.cpp`. Use this when you want to try some method using the `main.cpp` file. Note that we will not grade your `main.cpp` use it just for your own personal tests. To make and run use the following commands:
  ```
  make
  ./poly
  ```
- Making the test program. This will compile `polynomial.cpp` and `unit-test.cpp`. Use this to run the test script for this assignment. This will be the way that we will grade your lab. What you get on this test is going to be your raw-grade, that is your grade before the deductions that you may get for not following instructions, bad programming practices, etc. To make and run the test use the following commands:
  ```
  make test
  ./test
  ```

### Testing for Memory Leaks
To test for memory leaks we will use the test script. The following commands can be use:
```
make test
valgrind --leak-check=full ./test
```
Make sure that your code reports 0 errors, having memory leaks will affect your grade.

### Cleaning your Directory
If you want, you can run the command below to delete the binary code files (.o files and executables):
```
make clean
```

### Sample of the Command Line Output
```
root@84306837df10:/development/CSC2431/polynomial-solution# ll
total 36
drwxr-xr-x 13 root root  416 Apr  7 09:22 ./
drwxr-xr-x  6 root root  192 Apr  5 15:24 ../
drwxr-xr-x 12 root root  384 Apr  7 09:08 .git/
-rw-r--r--  1 root root  318 Apr  5 15:24 .gitignore
drwxr-xr-x  9 root root  288 Apr  7 09:09 .idea/
-rw-r--r--  1 root root  140 Apr  7 09:09 CMakeLists.txt
-rw-r--r--  1 root root 4198 Apr  5 15:24 README.md
drwxr-xr-x  8 root root  256 Apr  7 09:08 cmake-build-debug/
-rw-r--r--  1 root root  558 Apr  7 09:22 main.cpp
-rw-r--r--  1 root root  417 Apr  5 15:24 makefile
-rw-r--r--  1 root root 3333 Apr  5 15:28 polynomial.cpp
-rw-r--r--  1 root root  827 Apr  5 15:24 polynomial.h
-rw-r--r--  1 root root 2168 Apr  5 15:26 unit-test.cpp

root@84306837df10:/development/CSC2431/polynomial-solution# make
g++ -c main.cpp -o main.o -std=c++14 -Wall -g
g++ -c polynomial.cpp -std=c++14 -Wall -g
g++ main.o polynomial.o -o poly -std=c++14 -Wall -g

root@84306837df10:/development/CSC2431/polynomial-solution# ./poly
Running the Main Program
-2.50x^3 +0.00x^2 +2.00x^1 -1.00
3 -1 2 0 -2.5


+2.50x^3 -0.00x^2 -2.00x^1 +1.00

root@84306837df10:/development/CSC2431/polynomial-solution# make test
g++ -c unit-test.cpp -o test.o -std=c++14 -Wall -g
g++ test.o polynomial.o -o test -std=c++14 -Wall -g
root@84306837df10:/development/CSC2431/polynomial-solution# ./test
	SUM PASSED
	SUBTRACT PASSED
	MULTIPLY PASSED
	DERIVED PASSED
	EVALUATE PASSED
TOTAL TESTS TO PASS = 5
UNIT TEST GRADE = 100
 =======================
Extra Credit : : : : : : : : :
	INTEGRATE FAILED
	DIVIDE FAILED
UNIT TEST FINAL GRADE = 100

root@84306837df10:/development/CSC2431/polynomial-solution# valgrind --leak-check=full ./test
==240== Memcheck, a memory error detector
==240== Copyright (C) 2002-2017, and GNU GPL'd, by Julian Seward et al.
==240== Using Valgrind-3.15.0 and LibVEX; rerun with -h for copyright info
==240== Command: ./test
==240==
	SUM PASSED
	SUBTRACT PASSED
	MULTIPLY PASSED
	DERIVED PASSED
	EVALUATE PASSED
TOTAL TESTS TO PASS = 5
UNIT TEST GRADE = 100
 =======================
Extra Credit : : : : : : : : :
	INTEGRATE FAILED
	DIVIDE FAILED
UNIT TEST FINAL GRADE = 100
==240==
==240== HEAP SUMMARY:
==240==     in use at exit: 0 bytes in 0 blocks
==240==   total heap usage: 77 allocs, 77 frees, 76,612 bytes allocated
==240==
==240== All heap blocks were freed -- no leaks are possible
==240==
==240== For lists of detected and suppressed errors, rerun with: -s
==240== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```



