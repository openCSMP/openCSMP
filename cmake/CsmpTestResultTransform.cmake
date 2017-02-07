find_program(XSLT_PROCESSOR NAMES xsltproc Xalan msxsl)
if(NOT XSLT_PROCESSOR)
  message(FATAL_ERROR "no xslt processor found")
else(NOT XSLT_PROCESSOR)
  message(STATUS "Using ${XSPT_PROCESSOR} as XSLT processor")
endif(NOT XSLT_PROCESSOR)

get_filename_component(_result ${XSLT_PROCESSOR} NAME_WE)
if(_result STREQUAL "xsltproc")
  set(XSLT_PROCESSOR_ARGS_PATTERN 
    "<XSLT_PROCESSOR> -o <XML_OUTPUT> <XSL_CATALOG> <XML_INPUT>")
elseif(_result STREQUAL "Xalan")
  set(XSLT_PROCESSOR_ARGS_PATTERN 
    "<XSLT_PROCESSOR> -o <XML_OUTPUT> <XML_INPUT> <XSL_CATALOG>")  
elseif(_result STREQUAL "msxsl")
  set(XSLT_PROCESSOR_ARGS_PATTERN 
    "<XSLT_PROCESSOR> -o <XML_OUTPUT> <XML_INPUT> <XSL_CATALOG>")  
endif(_result STREQUAL "xsltproc")

separate_arguments(XSLT_PROCESSOR_ARGS_PATTERN)

function(CTEST_GET_RESULTDIR _out_results_dir _build_dir)
  file(STRINGS "${_build_dir}/Testing/TAG" _tag_content)
  list(GET _tag_content 0 _stamp)
  message(STATUS "Setting ${_out_results_dir} to ${_build_dir}/Testing/${_stamp} in parent scope")
  set(${_out_results_dir} "${_build_dir}/Testing/${_stamp}" PARENT_SCOPE)
endfunction(CTEST_GET_RESULTDIR)


function(CTEST_GENERATE_RESULTINDEX out_index_xml results_dir)
  set(_index_xml "${results_dir}/index.xml")
  file(GLOB _result_files "${results_dir}/*.xml")
  list(REMOVE_ITEM _result_files "${_index_xml}")
  
  file(WRITE "${_index_xml}"
    "<?xml version=\"1.0\" encoding=\"utf8\"?>\n"
    "<Files>\n")
  foreach(item ${_result_files})
    file(RELATIVE_PATH _result "${results_dir}" "${item}")
    file(APPEND "${_index_xml}"
      "\t<File>${_result}</File>\n")                            
  endforeach(item)
  file(APPEND "${_index_xml}"
    "</Files>")
    
  set(${out_index_xml} "${_index_xml}" PARENT_SCOPE)
endfunction(CTEST_GENERATE_RESULTINDEX)


function(CTEST_XSL_TRANSFORM  xsl_catalog input_xml output_xml)
  string(REPLACE "<XML_OUTPUT>" "${output_xml}" _xslt_command "${XSLT_PROCESSOR_ARGS_PATTERN}")
  string(REPLACE "<XML_INPUT>" "${input_xml}" _xslt_command "${_xslt_command}")
  string(REPLACE "<XSL_CATALOG>" "${xsl_catalog}" _xslt_command "${_xslt_command}")
  string(REPLACE "<XSLT_PROCESSOR>" "${XSLT_PROCESSOR}" _xslt_command "${_xslt_command}")
  
  message(STATUS "XSLT command: ${_xslt_command}")
  execute_process(
    COMMAND ${_xslt_command}    
    RESULT_VARIABLE _failed
    ERROR_VARIABLE _error
  )
endfunction(CTEST_XSL_TRANSFORM)
