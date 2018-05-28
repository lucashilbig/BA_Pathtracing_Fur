#!/bin/bash
directory=$(dirname $0)

green='\033[92m'
yellow='\033[93m'
none='\033[0m\033[39m' # No Color
bold='\033[1m'
blue='\033[34m'
red='\033[31m'

clear_print_head() {
    clear
    printf "$bold$green%s\n----------------\n$none" "KIRK build setup"
}

list_options() {
    printf "$bold%s$none\n" "Select your IDE"
    printf "$blue%-5s$none %-30s\n" 1 "XCode"
    printf "$blue%-5s$none %-30s\n" 2 "Unix Makefiles"
    printf "$blue%-5s$none %-30s\n" 3 "Eclipse with CDT"
}

get_list_selection() {
    while [[ true ]]
    do
        printf "$bold$blue:$none " >&2
        local selection
        read selection
        if [[ $selection =~ ^[-0-9]+$ && $selection > 0 && $selection < 4 ]]
        then
            echo $selection
            break
        fi
        printf "$bold$red%s$none\n" "The input \"$selection\" is not valid." >&2
    done
}

make_folder() {
# $1 -> IDE choice
    local folder
    case $1 in
        1)
            folder="Xcode_RT"
            ;;
        2)
            folder="CLI_RT"
            ;;
        3)
            folder="Eclipse_RT"
            ;;
    esac
    mkdir -p $folder
    # Just to be safe, use absolute path
    echo $directory/$folder
}

make_cmake_flag() {
# $1 -> IDE choice
    case $1 in
        1)
            echo "Xcode"
            ;;
        2)
            echo "Unix Makefiles"
            ;;
        3)
            echo "Eclipse CDT4 - Unix Makefiles"
            ;;
    esac
}

check_doxygen()
{
    local doxy
    printf "$blue%s$none" "Do you want to enable doxygen? (Y/n): " >&2
    read doxy
    case $doxy in
    n|N) echo "" ;;
    *) echo "-DBUILD_DOCUMENTATION=ON" ;;
    esac
}

quit() {
    if [[ $1 == 0 ]]
    then
#clear_print_head
        printf "$bold$green%s$none\n" "Setup finished successfully."
    else
        if [[ "$2" != "" ]]
        then
            printf "$bold$red%s$none\n" "$2"
        else
            printf "$bold$red%s$none\n" "Setup failed."
        fi
    fi
    printf "Press [ENTER] to close..."
    read _
    cd $directory
    clear
    exit 0
}

##############################################################################
###
### Main
###
##############################################################################
clear_print_head

dep="$PWD/dependencies"

if [ ! -d $dep ]; then
    quit 1 "Could not find dependencies folder. Please place it in the project root path as \"dependencies\"."
fi

export CVK_DEPENDENCIES_OSX=$dep

# Quit setup if cmake is not found.
! command -v cmake > /dev/null && quit 1 "Cmake not installed."

echo
list_options
echo
selected_item=$(get_list_selection)
folder=$(make_folder $selected_item)
clear_print_head
echo
doxygen_flag=$(command -v doxygen > /dev/null && check_doxygen)
cmake_flag=$(make_cmake_flag $selected_item)

clear_print_head
printf "$blue%s$none\n" "Generate for $cmake_flag"
[[ "$doxygen_flag" != "" ]] && printf "$blue%s$none\n" "Doxygen enabled"
echo

cd $folder
CC=gcc CXX=g++ cmake -G"$cmake_flag" -DCMAKE_BUILD_TYPE=Release -lstdc++fs "$doxygen_flag" $directory/../src/
status=$?
if [ $status -ne 0 ] ; then quit 1 "Failed to setup CMake project."; fi
quit 0
