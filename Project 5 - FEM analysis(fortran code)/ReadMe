# How to Run the Program

This program performs truss and frame analysis. Follow the steps below to run the program:
##Installation Requirements
Please note that the installation of Intel® Fortran Compiler is necessary to run this program. Ensure that it is correctly installed on your system before attempting to execute the program.

## Input

1. Open the text file named `TRUSS_DYN_IMPORT.txt` for truss analysis or `txt.FRAME_DYN_IMPORT` for frame analysis located in the program project folder.
2. Enter the structure and load specifications in the following format:

### Part 1: Program Title
Enter the title of the program.

### Part 2: General Specifications
Enter the following details in order:
- Number of nodes
- Number of elements
- Number of boundary conditions
- Number of groups of elements with common characteristics
- Number of characteristics of each section (1- modulus of elasticity, 2- cross-sectional area, 3- moment of inertia, 4- density)
- Number of nodes in each element
- Number of degrees of freedom
- Number of dimensions
- Force components in each node
- Dynamic analysis time frame
- Time step
- Whether the load is sinusoidal or cosine (cosine = 1, sine = 2)
- Alpha coefficient value
- Delta coefficient value
- Alpha value for damping
- Beta value for damping (according to the formula 𝛽 + 𝛼 =   )

### Part 3: Element Category Specifications
Enter the following details in order:
- Number of the element category with common features
- Modulus of elasticity
- Cross section
- Moment of inertia
- Density

### Part 4: Element Specifications
Enter the following details:
- Element number
- Number of the first node of the element
- Number of the second node of the element
- Number of the category of elements with common features

### Part 5: Node Specifications
Enter the following details:
- Node number
- X-coordinates of the node
- Y-coordinates of the node

### Part 6: Constraint Specifications 
Enter the following details:
- Node number 
- Constraint status in degree one and its value 
- Constraint status in degree two and its value 
- Constraint status in degree three (if any) and its value 

### Part 7: Load Specifications 
Enter the following details:
- Node number 
- Frequency and amplitude of load in each degree 

## Running The Program

After entering all specifications, run the program. The output file for truss analysis will be saved in `TRUSS_DYN_EXPORT.txt` and for frame analysis in `FRAME_DYN_EXPORT.txt`.
