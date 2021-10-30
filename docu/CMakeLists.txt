# doxygen installation on debian-like: apt install doxygen graphviz

find_package(
    Doxygen
    COMPONENTS dot
    OPTIONAL_COMPONENTS mscgen dia
)

if (DOXYGEN_FOUND)
    # general:
    set(DOXYGEN_USE_MDFILE_AS_MAINPAGE "${CMAKE_SOURCE_DIR}/README.md")
    set(DOXYGEN_USE_MATHJAX YES) # better formula rendering
    set(DOXYGEN_JAVADOC_AUTOBRIEF "YES") # for short '///' docstring support
    # set(DOXYGEN_EXTRACT_ALL YES) # Note: generates nice sub-sectioning for README.md, but also generates docu for all elements, including .cpp files
    set(DOXYGEN_BUILTIN_STL_SUPPORT YES)
    set(DOXYGEN_EXTRACT_LOCAL_CLASSES NO)

    # nicer graphics
    set(DOXYGEN_HAVE_DOT YES) # default anyway
    set(DOXYGEN_DOT_IMAGE_FORMAT svg) # use vector graphics, not png
    set(DOXYGEN_DOT_TRANSPARENT YES) # especially good for dark mode (at least in firefox, currently not chrome)

    # using graphics from docu/images
    set(DOXYGEN_IMAGE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/images/")
    
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

    # excluding examples/*.py to avoid polluting documented namespace, skipp examples/README.md due to not working file links:
    set(DOXYGEN_EXCLUDE_PATTERNS "*/python/examples/*.py;*/test/*;magnumaf.pyx;*/magnumaf/examples/README.md;")

    # adding examples to docu:
    set(DOXYGEN_EXAMPLE_PATH "${CMAKE_SOURCE_DIR}/python/examples/;${CMAKE_SOURCE_DIR}/magnumaf/examples/")
    # for files to be added, use '## @example' in python example or '/// @example' for cpp in script, better though to mention at other place in src with e.g. '## @example sp4.py' or '\example sp4.cpp'.

    # use modern CSS from https://github.com/jothepro/doxygen-awesome-css :
    if(TRUE) # TOGGLE_DARKMODE
      # generate header.html and footer.html files with custom commands in WORKIND_DIRECTORY
      execute_process(COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/css/darkmode/generate_html_header_and_footer_inserting_custom_commands.sh
        ${CMAKE_CURRENT_SOURCE_DIR}/css/darkmode/insert_header.html
        ${CMAKE_CURRENT_SOURCE_DIR}/css/darkmode/insert_footer.html
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      )

      set(DOXYGEN_HTML_EXTRA_STYLESHEET "css/doxygen-awesome.css;css/doxygen-awesome-sidebar-only.css;css/darkmode/doxygen-awesome-sidebar-only-darkmode-toggle.css")
      set(DOXYGEN_GENERATE_TREEVIEW YES)
      set(DOXYGEN_HTML_EXTRA_FILES "css/darkmode/doxygen-awesome-darkmode-toggle.js")
      set(DOXYGEN_HTML_HEADER "${CMAKE_CURRENT_BINARY_DIR}/header.html")
      set(DOXYGEN_HTML_FOOTER "${CMAKE_CURRENT_BINARY_DIR}/footer.html")
    else()
      set(DOXYGEN_HTML_EXTRA_STYLESHEET "css/doxygen-awesome.css;css/doxygen-awesome-sidebar-only.css")
      set(DOXYGEN_GENERATE_TREEVIEW YES)
    endif() # TOGGLE_DARKMODE


    # make latex output and cmake target if found
    find_package(LATEX COMPONENTS PDFLATEX)
    if(LATEX_FOUND)
        set(DOXYGEN_LATEX_CMD_NAME "latex")
        set(DOXYGEN_GENERATE_LATEX YES)
        set(DOXYGEN_EXTRA_PACKAGES "{amsmath};{bm}")
        add_custom_target(
            docu_pdf
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
