Considering you are using a Linux-based OS, install git and the pre-requisites for g++ with

`sudo apt-get install git `

`sudo apt-get install build-essential` 

`sudo apt-get install cmake `

`sudo apt-get install autoconf`

Clone the library with the command

`git clone --recurse-submodules https://gitlab.com/palisade/palisade-release`

If the git command gives some errors, alternatively you can download the library from https://gitlab.com/palisade/palisade-release, go to the git pages of the third-party applications and manually download and extract the source codes to the related folder under **/third-party** (You can find the URLs in the file **.gitmodules**)

Create an empty '**build**' folder in the cloned palisade-release folder with

`mkdir build`

Then go inside and run cmake as

`cd build`

`cmake ..`

After that, just run the commands below, and this should conclude the installation.

`make`

`sudo make install  `

More detailed instructions to install the PALISADE can be found under the “Getting Started with PALISADE” section in https://gitlab.com/palisade/palisade-release/-/wikis/home.
