# imagemicro
Small tkinter gui tool for basic greyscale image processing using PIL and scipy
<img src="thresholddemo.jpg" style="width: 500px;"/>

Usage: A file may be opened using the file dialog. This is placed in the primary greyscale buffer. Thresholding places the result of the thresholding operation in the primary binary buffer. This is superimposed on the primary greyscale buffer in red. The primary greyscale buffer may be saved into one of the greyscale slots for later use. Similarly, the primary binary buffer may be saved into one of the binary slots for later use. GS Operations (greyscale operations) include some routines for edge detection and the watershed algorithm. Bin Operations (binary operations) operate on binary images. The binary image may be inverted, or a logical operation (AND, OR, XOR) may be used. This performs the chosen operation on two of the binary slots (default: 1 and 2) and places the result into the primary binary buffer. From scipy.ndimage come binary operations erosion, dilation, opening, closing and fill holes. The first four open a menu which permit setting the size of the kernel. Currently, this is a a circle kernel. A logger is included which prints the python code used to produce the actions performed in the gui by the user.

Todo: Clean up the code, separate UI from function further. Features that would be useful are manual removal of binary objects, undo functionality, binary masking of greyscale images and greyscale image arithmetic. Cool features would be low/high pass filters.


Please note that this is a work in progress and primarily created as a learning exercise.

## Usage
The easiest way to use the program is through virtualenv. Make sure you have python (2.7) and virtualenv installed. Clone the repository and, in the repository folder, run the following commands.

```bash
virtualenv env
source ./env/bin/activate
pip2 install -r requirements.txt
```

Now, run `python2 main.py` to start the program.
