/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "StringBasics.h"


int main(int argc, char ** argv)
{
    String outputString = "Hello, welcome to SampleProgram.\n";

    IFILE filePtr;
    
    filePtr = ifopen("-", "w");

    if(filePtr == NULL)
    {
        std::cerr << "Failed to open stdout\n";
        return(-1);
    }

    ifwrite(filePtr, outputString.c_str(), outputString.Length());

    ifclose(filePtr);
    return(0);
}