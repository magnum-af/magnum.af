# doxygen installation on debian-like: apt install doxygen graphviz

find_package(
    Doxygen
    REQUIRED dot
    OPTIONAL_COMPONENTS mscgen dia
)

if (DOXYGEN_FOUND)
    # general:
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_SOURCE_DIR}/README.md")
    set(DOXYGEN_USE_MATHJAX YES) # better formula rendering
    set(DOXYGEN_JAVADOC_AUTOBRIEF "YES") # for short '///' docstring support
    set(DOXYGEN_EXTRACT_LOCAL_CLASSES "NO")
    
    # resetting defaults to make sure:
    set(DOXYGEN_RECURSIVE "YES") # is default anyway
    set(DOXYGEN_TOC_INCLUDE_HEADINGS "5") # is default anyway

    # make file links relative:
    set(DOXYGEN_FULL_PATH_NAMES "YES") # Note: is default anyway
    set(DOXYGEN_STRIP_FROM_PATH "${CMAKE_SOURCE_DIR}")

    # document .pyx as python
    # Note: doxygen would mix c++ and python namespace 'magnumaf' and all classes within,
    # as python namespace is infered from filename. (c++: magnumaf::, py: magnumaf.py).
    # Therefore we add a symlink called magnumafpy.pyx and exclude the original magnumaf.pyx:
    # TODO remove symlink, rename file with cmake only in build dir!
    set(DOXYGEN_EXTENSION_MAPPING ".pyx=python")
    set(DOXYGEN_FILE_PATTERNS *.cpp *.hpp *.md *.py *.pyx) # Note: don't use quotes here
    set(DOXYGEN_EXCLUDE_PATTERNS "*/test/*;*/examples/*;magnumaf.pyx")

    # make latex output and cmake target if found
    find_package(LATEX COMPONENTS PDFLATEX)
    if(LATEX_FOUND)
        set(DOXYGEN_LATEX_CMD_NAME "latex")
        set(DOXYGEN_GENERATE_LATEX YES)
        set(DOXYGEN_EXTRA_PACKAGES "{amsmath};{bm}")
        add_custom_target(
            docu_latex
            COMMAND make
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/latex"
            COMMENT "Building doxygen latex pdf documentation."
        )
    endif(LATEX_FOUND)

    # generate Doxyfile with above settings
    doxygen_add_docs(
        docu_html ${CMAKE_SOURCE_DIR}
        COMMENT "Generating doxygen html documentation."
    )

endif(DOXYGEN_FOUND)
