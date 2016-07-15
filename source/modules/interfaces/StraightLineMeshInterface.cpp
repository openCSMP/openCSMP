#include "StraightLineMeshInterface.h"
#include "Mesher1D.h"
#include "TextFileInterface.h"

using namespace std;

namespace csmp {

template<size_t dim>
StraightLineMeshInterface<dim>
::StraightLineMeshInterface( )
{
}

template<size_t dim>
StraightLineMeshInterface<dim>
::StraightLineMeshInterface( csmp::VSet<dim>& vset,
                             const std::string& mesh_file )
{
     Initialize( vset, mesh_file );
}

// MODES

template<size_t dim>
void StraightLineMeshInterface<dim>
::WriteInfoFile()
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    if( error_handler.Verbose() )
        cout <<"\n\nStraightLineMeshInterface::WriteInfoFile: Write an info file"<< endl;

    ofstream  ofs;
    std::string file_extension("-info.txt");

    // 0. Open file
    csmp::openFile( ofs, "StraightLineMeshInterface", file_extension );

    // 1. Write info file
    this->WriteSampleFile( ofs );

    ofs.close();

    if( error_handler.Verbose() )
        cout <<"\nStraightLineMeshInterface::WriteInfoFile: Writing of info file is completed."<< endl;

}

template<size_t dim>
void StraightLineMeshInterface<dim>
::WriteSampleFile( std::ofstream& ofs)
{
    ofs <<   "## " << "====================================================================";
    ofs << "\n## " << "Sample mesh file";

    ofs << "\n";

    ofs << "\n####### " << "IMPORTANT NOTE 1: spaces, tabs and other delimiters";
    ofs << "\n## " << "Tabs, commas or semicolons are used as column data separators, but";
    ofs << "\n## " << "DO NOT USE SPACES AS SEPARATORS!!!!";
    ofs << "\n## " << "Spaces can only be used in object names and not as separators.";
    ofs << "\n## " << "Using spaces in inpropper way can cause unpredicted problems.";
    ofs << "\n## " << "Please double check your configuration file before using it.";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n## " << "name of object\tdata1\ttdata 2 as some name of smth.\tdata3";

    ofs << "\n";

    ofs << "\n####### " << "IMPORTANT NOTE 2: comments";
    ofs << "\n## " << "A comment line should be preceded by the sharp sign '#'";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n## " << "# comment line";

    ofs << "\n";

    ofs << "\n## " << "Please, find bellow the description of the main keywords";
    ofs << "\n## " << "which you can use for configuration of your application.";
    ofs << "\n## " << "====================================================================";

    ofs << "\n";

    ofs << "\n## " << "====================================================================";
    ofs << "\n## " << "Keywords from StraightLineMeshInterface class";
    ofs << "\n## " << "====================================================================";

    ofs << "\n";

    ofs << "\n####### " << "KEYWORD: *** UNIFORM_MESH ***";
    ofs << "\n####### " << "DESCRIPTION:";
    ofs << "\n## " << "Allows to build uniform mesh.";
    ofs << "\n## " << "One need to put the length of the domain and number of desired elements.";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n## " << "UNIFORM_MESH\tlength\tnumber of elements";

    ofs << "\n";

    ofs << "\n####### " << "KEYWORD: *** REFINED_MESH ***";
    ofs << "\n####### " << "DESCRIPTION:";
    ofs << "\n## " << "Allows to build refined mesh.";
    ofs << "\n## " << "One need to put the length of the domain,";
    ofs << "\n## " << "minimum and maximum allowed sizes of elements,";
    ofs << "\n## " << "width of transition zone, where mesh size is gradually increases, and";
    ofs << "\n## " << "density function, which by default is linear.";
    ofs << "\n## " << "Avaliable density functions are:";
    ofs << "\n## " << "lin: a * x + b, defaults: a = 1.0, b = 0.0";
    ofs << "\n## " << "exp: b * exp( a * x ), defaults: a = 1.0, b = 1.0";
    ofs << "\n## " << "erf: b * erf( a * x ), defaults: a = 3.0, b = 1.0";
    ofs << "\n## " << "The error function's table is: erf(0.0)=0.0, erf(0.5)~0.5, erf(1.0)~0.8, erf(3.0)~1.0";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n# " << "# refined mesh with density function: x";
    ofs << "\n## " << "REFINED_MESH\tlength\tdx_min\tdx_max\twidth of transition zone";
    ofs << "\n# " << "# refined mesh with denisty function: 3*x+1";
    ofs << "\n## " << "REFINED_MESH\tlength\tdx_min\tdx_max\twidth of transition zone\tlin\t3.0\t1.0";
    ofs << "\n# " << "# refined mesh with density function: exp(x)";
    ofs << "\n## " << "REFINED_MESH\tlength\tdx_min\tdx_max\twidth of transition zone\texp";
    ofs << "\n# " << "# refined mesh with density function: 2*exp(5*x)";
    ofs << "\n## " << "REFINED_MESH\tlength\tdx_min\tdx_max\twidth of transition zone\texp\t5.0\t2.0";
    ofs << "\n# " << "# refined mesh with density function: erf(3*x)";
    ofs << "\n## " << "REFINED_MESH\tlength\tdx_min\tdx_max\twidth of transition zone\terf";
    ofs << "\n# " << "# refined mesh with density function: 3*erf(2*x)";
    ofs << "\n## " << "REFINED_MESH\tlength\tdx_min\tdx_max\twidth of transition zone\terf\t2.0\t3.0";

    ofs << "\n";

    ofs << "\n####### " << "KEYWORD: *** CUSTOM_MESH ***";
    ofs << "\n####### " << "DESCRIPTION:";
    ofs << "\n## " << " Allows to build mesh based on provided coordinates of nodes.";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n## " << "CUSTOM_MESH";
    ofs << "\n## " << "node1_x\tnode1_y\tnode1_z";
    ofs << "\n## " << "node2_x\tnode2_y\tnode2_z";
    ofs << "\n## " << "node3_x\tnode3_y\tnode3_z";
    ofs << "\n## " << "node4_x\tnode4_y\tnode4_z";
    ofs << "\n## " << "node5_x\tnode5_y\tnode5_z";

    ofs << "\n";

    ofs << "\n####### " << "KEYWORD: *** CORNER_POINTS ***";
    ofs << "\n####### " << "DESCRIPTION:";
    ofs << "\n## " << "Allows to adjust the existing line mesh to the provided corner points.";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n## " << "CORNER_POINTS";
    ofs << "\n## " << "corner1_x\tcorner1_y\tcorner1_z";
    ofs << "\n## " << "corner2_x\tcorner2_y\tcorner2_z";

    ofs << "\n";

    ofs << "\n####### " << "KEYWORD: *** SPLITNODES ***";
    ofs << "\n####### " << "DESCRIPTION:";
    ofs << "\n## " << "Allows to add splitted nodes at provided locations.";
    ofs << "\n####### " << "USAGE:";
    ofs << "\n## " << "SPLITNODES";
    ofs << "\n## " << "point1_x\tpoint1_y\tpoint1_z";
    ofs << "\n## " << "point2_x\tpoint2_y\tpoint2_z";
    ofs << "\n## " << "point3_x\tpoint3_y\tpoint3_z";
}


/** Builds Model based on the information in text file
  @author R. Manasipov
*/

template<size_t dim>
void StraightLineMeshInterface<dim>
::Initialize( csmp::VSet<dim>& vset, const std::string& fname )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    // Identify pointer to vset
    vset_ = &vset;

    if( error_handler.Verbose() )
        cout <<"\nStraightLineMeshInterface<"<< dim <<">::Initialize: Build model from file..."<< endl;

    ifstream  ifs;
    string    file_header;
    size_t    line_length(256);

    // 0. Open the file
    std::string filename = fname;
    filename += "-mesh.txt";
    csmp::openFile( ifs, filename );

    // 1. Read the title
    csmp::readFileHeader( ifs, file_header, error_handler.Verbose() );

    // 2. Read the data
    csmp::TextFileInterface<StraightLineMeshInterface<dim> > reader( line_length );

    reader.AddKeyword("UNIFORM_MESH" );
    reader.AddKeyword("REFINED_MESH" );
    reader.AddKeyword("CUSTOM_MESH"  );
    reader.AddKeyword("CORNER_POINTS");
    reader.AddKeyword("SPLITNODES"   );

    reader.ReadFile( ifs, this, &csmp::isCommentLine, (&StraightLineMeshInterface<dim>::ReadMesh) );

    if( error_handler.Verbose() )
        cout <<"\n\nStraightLineMeshInterface<"<< dim <<">::Initialize: Mesh is Builded!."<< endl;
}

/// Generic function which reads blocks of data marked by keyword

template<size_t dim>
bool StraightLineMeshInterface<dim>
::ReadMesh( std::ifstream& ifs, char* text_line, size_t line_length,
            std::set<std::string>& keywords, std::string& keyword, std::vector<std::string>& keyword_params )
{
    csmp::ErrorHandler& error_handler( csmp::ErrorHandler::Instance() );

    csmp::Mesher1D<dim> mesher;
    if ( keyword == "UNIFORM_MESH" )
    {
        if( keyword_params.size() < 2 )
            return false;
        const double  length  ( atof(keyword_params[0].c_str() ) );
        const size_t  elements( atoi(keyword_params[1].c_str() ) );
        mesher.BuildUniformMesh( *vset_, length, elements );
    }
    else if( keyword == "REFINED_MESH" )
    {
        if( keyword_params.size() < 4 )
            return false;
        double length( atof(keyword_params[0].c_str() ) );
        double dx_min( atof(keyword_params[1].c_str() ) );
        double dx_max( atof(keyword_params[2].c_str() ) );
        double width_of_transition_zone( atof(keyword_params[3].c_str() ) );
        csmp::LinDensity   lindens;
        csmp::ExpDensity   expdens;
        csmp::ErfDensity   erfdens;
        csmp::MeshDensity* meshdens = static_cast<csmp::MeshDensity*>(&lindens);
        if( keyword_params.size() > 4 )
        {
            if( keyword_params[4] == std::string("lin") )
            {
                if( keyword_params.size() > 5 )
                    lindens.SetA( atof(keyword_params[5].c_str() ) );
                if( keyword_params.size() > 6 )
                    lindens.SetB( atof(keyword_params[6].c_str() ) );
                meshdens = static_cast<csmp::MeshDensity*>(&lindens);
            }
            else if( keyword_params[4] == std::string("exp") )
            {
                if( keyword_params.size() > 5 )
                    expdens.SetA( atof(keyword_params[5].c_str() ) );
                if( keyword_params.size() > 6 )
                    expdens.SetB( atof(keyword_params[6].c_str() ) );
                meshdens = static_cast<csmp::MeshDensity*>(&expdens);
            }
            else if( keyword_params[4] == std::string("erf") )
            {
                if( keyword_params.size() > 5 )
                    erfdens.SetA( atof(keyword_params[5].c_str() ) );
                if( keyword_params.size() > 6 )
                    erfdens.SetB( atof(keyword_params[6].c_str() ) );
                meshdens = static_cast<csmp::MeshDensity*>(&erfdens);
            }
        }
        if( width_of_transition_zone == 0.0 )
            mesher.BuildUniformMesh( *vset_, length, static_cast<size_t>(length/dx_min) );
        else
            mesher.BuildRefinedMesh( *vset_, length, dx_min, dx_max, width_of_transition_zone, meshdens );
    }
    else if( keyword == "CUSTOM_MESH" )
    {
        char*       token(0);
        const char* delims =" :,\t,\n,\r";
        std::set<csmp::Point<dim> > nodes;
        csmp::Point<dim> p( 0.0 );
        do {
            if ( !csmp::isCommentLine(text_line) )
            {
                // read x coordinate
                token = strtok( text_line, delims );
                p[0] = atof( token );
                // read remaining coordinates or assign it to 0.0 if notning is present
                for ( size_t dimension=1U; dimension < dim; dimension++ )
                {
                    token = strtok( NULL, delims );
                    if ( token != NULL )
                        p[ dimension ] = atof(token);
                    else
                    {
                        for ( size_t i = dimension; i<dim; i++ )
                            p[ i ] = 0.0;
                        break;
                    }
                }
                nodes.insert( p );
            }
            ifs.getline( text_line, line_length );
        }
        while ( !csmp::isBlankLine(text_line) && !ifs.eof() );

        std::vector<csmp::Point<dim> > node_vector( nodes.begin(), nodes.end() );
        mesher.BuildCustomMesh( *vset_, node_vector );
    }
    else if( keyword == "CORNER_POINTS" )
    {
        char*       token(0);
        const char* delims =" :,\t,\n,\r";
        std::set<csmp::Point<dim> > nodes;
        csmp::Point<dim> p( 0.0 );
        size_t counter( 0 );
        do {
            if ( !csmp::isCommentLine(text_line) )
            {
                if( counter < 2 )
                {
                    // read x coordinate
                    token = strtok( text_line, delims );
                    p[0] = atof( token );
                    // read remaining coordinates or assign it to 0.0 if notning more is provided
                    for ( size_t dimension=1U; dimension < dim; dimension++ )
                    {
                        token = strtok( NULL, delims );
                        if ( token != NULL )
                            p[ dimension ] = atof(token);
                        else
                        {
                            for ( size_t i = dimension; i < dim; i++ )
                                p[ i ] = 0.0;
                            break;
                        }
                    }
                    nodes.insert( p );
                    counter++;
                }
            }
            ifs.getline( text_line, line_length );
        }
        while ( !csmp::isBlankLine(text_line) && !ifs.eof() );

        if( counter != 2 )
            return false;

        mesher.AssignCornerPoints( *vset_, *nodes.begin(), *nodes.rbegin() );
    }
    else if( keyword == "SPLITNODES" )
    {
        char*       token(0);
        const char* delims =" :,\t,\n,\r";
        std::set<csmp::Point<dim> > nodes;
        csmp::Point<dim> p( 0.0 );
        do {
            if ( !csmp::isCommentLine(text_line) )
            {
                // read x coordinate
                token = strtok( text_line, delims );
                p[0] = atof( token );
                // read remaining coordinates or assign it to 0.0 if notning more is provided
                for ( size_t dimension=1U; dimension < dim; dimension++ )
                {
                    token = strtok( NULL, delims );
                    if ( token != NULL )
                        p[ dimension ] = atof(token);
                    else
                    {
                        for ( size_t i = dimension; i<dim; i++ )
                            p[ i ] = 0.0;
                        break;
                    }
                }
                nodes.insert( p );
            }
            ifs.getline( text_line, line_length );
        }
        while ( !csmp::isBlankLine(text_line) && !ifs.eof() );

        mesher.InsertSplitNodes( *vset_, nodes );
    }
    else
    {
        std::string msg = "Unknown keyword: ";
        msg += keyword;
        error_handler.notice( csmp::FATAL_ERROR,
                              "StraightLineMeshInterface<dim>::ReadMesh",
                              msg
                            );
        return false;
    }

    return true;
}




template class StraightLineMeshInterface<1U>;
template class StraightLineMeshInterface<2U>;
template class StraightLineMeshInterface<3U>;

} // end namespace csmp
