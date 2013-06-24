// /********************************************\
//  * Author: Marc Khoury
//  * Application: IsoRender
//  * An application for viewing isosurface meshes
// \********************************************/

/*
Copyright (c) 2010 Marc Khoury

Please send bug reports to marc.khry@gmail.com

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
*/

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <memory.h>
#include <GL/glut.h>
#include <GL/glu.h>
#include <GL/gl.h>
#include <GL/glui.h>
#include "Vector.h"
#include "ijkmcubeIO.h"

using namespace IJK;
using namespace IJKMCUBE;
using namespace std;

#define ESC 27
#define MAX(a,b) ((a) > (b) ? (a) : (b))

//Initialize Opengl
void init();

//GLUT callback functions
void display();
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void reshape(int w, int h);

//User Defined Operations
void set_ambient();
void set_light_positions();
void calc_world_coords(double &px,double &py,double &pz, int x, int y, int *viewport);
void setup_glui();
void setup_marching_cubes(int id);
void create_isosurface();
void process_isosurface();
void process_normals(COORD_ARRAY& coord, VERTEX_INDEX_ARRAY& simplices, int dimension);
void parse_command_line(int argc, char* argv[]);
void read_data();
void clip(int id);
void center_on_screen();
void memory_exhaustion();
void draw_bounding_box();
void reset(int id);
void write_off(int id);
void write_obj(int id);
void write_ppm(int id);
void construct_isosurface(const IO_INFO & info, const MC_DATA & mc_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);
void construct_interval_volume(const IO_INFO & info, const MC_DATA & mc_data, MCUBE_TIME & mcube_time, IO_TIME & io_time);
void get_matrix();
void inverse(const GLdouble *m, GLdouble *out );

//Global variables
int view_angle = 60;
float light_position[4] = {0,0,-1,1}; //Position of positional light source
int mouse_x = 0 ,mouse_y = 0;
double drag_x = 0, drag_y = 0, drag_z = 0;
bool mouse_left = false, mouse_right = false, mouse_middle = false; //Current mouse button being pressed
double w_left = 0.0, w_right = 0.0, w_bottom = 0.0, w_top = 0.0, z_near = 0.1, z_far = 100; //Used for orthographic projection into world coordinates
float bounding_box_vertices[8][3] = {{1,1,1},{-1,1,1},{-1,-1,1},{1,-1,1},{1,1,-1},{-1,1,-1},{-1,-1,-1},{1,-1,-1}}; //Bounding box defined:(-1,-1,-1)(1,1,1)
float bounding_box_faces[6][4] = {{0,1,2,3},{7,6,5,4},{4,0,3,7},{3,2,6,7},{1,5,6,2},{0,4,5,1}}; //Faces for bounding box

//Variables specific to GLUI and marching cubes set up
int main_window;
int glui_window_id;
int window_width = 930;
int window_height = 600;
int window_x, window_y;
int front_face_index;
int back_face_index;
int dim; //Dimension of the scalar grid
int subsample, supersample; //subsample and supersample resolution
float isovalue = 0.0;
float clip_plane = 0.1;
int draw_mode; //drawing primitive
int draw_bb = 0; //Draw bounding box
int back_face = 0;
bool recalculate = true; //Set to true when new isosurface is generated
char window_title[] = "IsoRender";
GLUI *glui_window;
Vector *normals; //Array of normal vectors
Point min_p, max_p; //Min and max points for the bounding box
MC_ISOSURFACE isosurface; //Current isosurface
MC_SCALAR_GRID full_scalar_grid; //representation of scalar grid
MC_DATA mc_data; //mc data structures and flags
NRRD_INFO nrrd_info; //Data read from nrrd input file
IO_INFO io_info; //Settings for Marching Cubes algorithm
GLUI_Spinner *isovalue_spinner;
//Matrices
double matrix1[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};          
double matrix2[16] = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};  
float front_colors[7][4] = {{0,1,0,1},{1,0,0,1},{0,0,1,1},{1,1,0,1},{0,1,1,1},{1,0,1,1},{1,1,1,1}};
float back_colors[7][4] = {{1,0,1,1},{1,1,0,1},{0,1,1,1},{0,1,0,1},{1,0,0,1},{0,0,1,1},{1,1,1,1}};

int main(int argc, char* argv[])
{
	parse_command_line(argc, argv);
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB  |GLUT_DEPTH);
	glutInitWindowSize(window_width, window_height);
	center_on_screen();
	glutInitWindowPosition(window_x,window_y);
	main_window = glutCreateWindow("IsoRender");
	init();

	glutDisplayFunc(display);  //Set glut call backs, I should probably conform to glui callbacks but this seems to work fine
	glutKeyboardFunc(keyboard);
	glutMotionFunc(motion);
	glutMouseFunc(mouse);
	glutReshapeFunc(reshape);

	setup_glui(); //Set up glui
	glutMainLoop();
	return 0;
}

void init()
{
	glClearColor(0.0,0.71,1.0,1.0);
	glShadeModel(GL_SMOOTH);
	
	float ambient[4] = {0,0,0,1};
	float diffuse[4] = {1,1,1,1};
	float specular[4] = {1,1,1,1};
	glEnable(GL_LIGHTING); //Set up lighting parameters
	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); //Enable two sided lighting so we can color back faces different colors
	glEnable(GL_DEPTH_TEST); //Enable depth test
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(view_angle, 1, z_near, z_far); //Define view frustrum
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void clip(int id)
{
	glViewport(0,0,window_width,window_height);
	w_top    =  1.0;
	w_bottom = -1.0;
	w_left   = -(double)window_width/(double)window_height;
	w_right  = -w_left;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
    	gluPerspective(view_angle,(double)window_width/(double)window_height, clip_plane, z_far); //Set the new aspect ratio

	glMatrixMode(GL_MODELVIEW);
	z_near = clip_plane;
}

void setup_glui()
{
	glui_window = GLUI_Master.create_glui_subwindow(main_window,GLUI_SUBWINDOW_BOTTOM); //Creat a GLUI subwindow to the left of the glut viewing window
	glui_window_id = glui_window->get_glut_window_id(); //Get glut window id
	glui_window->set_main_gfx_window(main_window);
	GLUI_Panel *left_panel = glui_window->add_panel(window_title); //Create a raised panel to place widgets
	isovalue_spinner = glui_window->add_spinner_to_panel(left_panel, "Isovalue", GLUI_SPINNER_FLOAT, &isovalue, 0, setup_marching_cubes);
	isovalue_spinner->set_float_limits(0.0, 1000.0); //Set limit on isovalue [0,1000]
	GLUI_Spinner *clipping_plane = glui_window->add_spinner_to_panel(left_panel, "Clip", GLUI_SPINNER_FLOAT, &clip_plane, 0, clip);
	clipping_plane->set_float_limits(0.1, 10.0); //Set limit on clipping plane	
	glui_window->add_column_to_panel(left_panel, true);		
	GLUI_Spinner *subsample_spinner = glui_window->add_spinner_to_panel(left_panel, "Subsample", GLUI_SPINNER_INT, &subsample, 0, setup_marching_cubes);
	subsample_spinner->set_int_limits(1, 50); //Subsample spinner
	GLUI_Spinner *supersample_spinner = glui_window->add_spinner_to_panel(left_panel, "Supersample", GLUI_SPINNER_INT, &supersample, 0, setup_marching_cubes);
	supersample_spinner->set_int_limits(1, 50); //Supersample spinner
	glui_window->add_column_to_panel(left_panel, true);
	GLUI_Listbox *front_color_list= glui_window->add_listbox_to_panel(left_panel,"Front ", &front_face_index);
	front_color_list->add_item(0,"Green");
	front_color_list->add_item(1,"Red");
	front_color_list->add_item(2,"Blue");
	front_color_list->add_item(3,"Yellow");
	front_color_list->add_item(4,"Cyan");
	front_color_list->add_item(5,"Purple");
	front_color_list->add_item(6,"White");
	GLUI_Listbox *back_color_list = glui_window->add_listbox_to_panel(left_panel,"Back ",&back_face_index);
	back_color_list->add_item(0,"Purple");	
	back_color_list->add_item(1,"Yellow");
	back_color_list->add_item(2,"Cyan");
	back_color_list->add_item(3,"Green");
	back_color_list->add_item(4,"Red");
	back_color_list->add_item(5,"Blue");	
	back_color_list->add_item(6,"White");
	glui_window->add_column_to_panel(left_panel, true);
	GLUI_RadioGroup * draw_mode_radio = glui_window->add_radiogroup_to_panel(left_panel, &draw_mode); //Radio group for drawing solid or wireframe
	glui_window->add_radiobutton_to_group(draw_mode_radio, "Solid Surface");
	glui_window->add_radiobutton_to_group(draw_mode_radio, "Wire Frame");
	glui_window->add_column_to_panel(left_panel, true);
	glui_window->add_checkbox_to_panel(left_panel,"Draw Bounding Box", &draw_bb); //checkbox for bounding box
	glui_window->add_checkbox_to_panel(left_panel,"Color Back Faces", &back_face); //checkbox for back face seperate color
	glui_window->add_column_to_panel(left_panel, true);
	glui_window->add_button_to_panel(left_panel,"Write OFF File",1,write_off); //Write OFF file callback
	glui_window->add_button_to_panel(left_panel,"Write OBJ File",1,write_obj); //Write OBJ file callback
	glui_window->add_column_to_panel(left_panel, true);
	glui_window->add_button_to_panel(left_panel,"Write PPM File",1,write_ppm); //write ppm file, still buggy	
	glui_window->add_button_to_panel(left_panel,"Reset",1,reset); //Set up reset callback and button on panel
}

void display()
{
	GLUI_Master.auto_set_viewport();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	set_ambient();
	set_light_positions();
	glEnable(GL_LIGHT0); //Enable lights
	if(recalculate)
	{process_isosurface();} //Need to transform mesh and calculate normals
	if(draw_mode == 1) //Set wireframe mode
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
	if(back_face == 0) //if we want to color the back faces a different color from the front faces
	{
		float zero[4] = {0,0,0,1};
		float specular[4] = {0,0,0,1};
		float shininess[1] = {1};
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, front_colors[front_face_index]);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, zero);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shininess);
	}
	else
	{	
		float zero[4] = {0,0,0,1};
		float specular[4] = {0,0,0,1};
		float shininess[1] = {1};
		glMaterialfv(GL_FRONT, GL_DIFFUSE, front_colors[front_face_index]);
		glMaterialfv(GL_FRONT, GL_AMBIENT, zero);
		glMaterialfv(GL_FRONT, GL_SPECULAR, specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, shininess);
		glMaterialfv(GL_BACK, GL_DIFFUSE, back_colors[back_face_index]);
		glMaterialfv(GL_BACK, GL_AMBIENT, zero);
		glMaterialfv(GL_BACK, GL_SPECULAR, specular);
		glMaterialfv(GL_BACK, GL_SHININESS, shininess);
	}
	glPushMatrix();
	glLoadIdentity();
	glTranslatef(0,0,-4); //Tranlate away from default location of camera
	glMultMatrixd(matrix1);	 //Apply transformation determined by user input
	if(draw_bb == 1)
	{draw_bounding_box();}
	for(int js = 0; js < isosurface.simplex_vert.size(); js = js + dim) //Render the mesh with proper normals
	{
		glBegin(GL_TRIANGLES);
		for(int k = 0; k < dim; k++)
		{
			int iv = isosurface.simplex_vert[js + k];
			glNormal3f(normals[iv][0],normals[iv][1],normals[iv][2]);			
			glVertex3f(isosurface.vertex_coord[iv*dim],isosurface.vertex_coord[iv*dim+1], isosurface.vertex_coord[iv*dim+2]);
		}
		glEnd();
	}
	glTranslatef(0,0,4); //Translate back to origin
	glPopMatrix(); 
	glutSwapBuffers();
}

void process_isosurface()
{
	if(isosurface.vertex_coord.size() > 0) //Set inital max and min values
	{
		max_p.x = isosurface.vertex_coord[0];
		max_p.y = isosurface.vertex_coord[1];
		max_p.z = isosurface.vertex_coord[2];
		min_p.x = isosurface.vertex_coord[0];
		min_p.y = isosurface.vertex_coord[1];
		min_p.z = isosurface.vertex_coord[2];
	}
	for(int js = dim; js < isosurface.vertex_coord.size(); js = js + dim)  //Find the maximum and minimum points on the surface
	{
		float x = isosurface.vertex_coord[js];
		float y = isosurface.vertex_coord[js+1];
		float z = isosurface.vertex_coord[js+2];
		min_p.x = x < min_p.x ? x : min_p.x;
		min_p.y = y < min_p.y ? y : min_p.y;
		min_p.z = z < min_p.z ? z : min_p.z;
		max_p.x = x > max_p.x ? x : max_p.x;
		max_p.y = y > max_p.y ? y : max_p.y;
		max_p.z = z > max_p.z ? z : max_p.z;
	}
	float center[3] = {(max_p.x+min_p.x)/2,(max_p.y+min_p.y)/2,(max_p.z+min_p.z)/2}; //Calculate the center of the bounding box
	float scale =  2.0f/(MAX((max_p.x-min_p.x), MAX((max_p.y-min_p.y),(max_p.z-min_p.z)))); //Calculate scale needed to place mesh in unit cube
	for(int js = 0; js < isosurface.vertex_coord.size(); js = js + dim) //Translate all vertices to the origin and scale them to fit in cube
	{
		isosurface.vertex_coord[js] -= center[0];
		isosurface.vertex_coord[js+1] -= center[1];
		isosurface.vertex_coord[js+2] -= center[2];
		isosurface.vertex_coord[js] *= scale;
		isosurface.vertex_coord[js+1] *= scale;
		isosurface.vertex_coord[js+2] *= scale;	
	}

	process_normals(isosurface.vertex_coord,isosurface.simplex_vert, dim); //Calculate new normals for the resulting surface
}


//------------------------------
//Draw the bounding box defined
//by the unit cube
//------------------------------
void draw_bounding_box()
{
	glDisable(GL_LIGHTING);  //Turn off lighting so that the bounding box will always be a white color visible from every angle
	glColor3f(1,1,1);
	for(int i = 0; i < 6; i++) 
	{
		glBegin(GL_LINE_LOOP);
		for(int j = 0; j < 4; j++)
		{
			int vertex = bounding_box_faces[i][j];
			glVertex3f(bounding_box_vertices[vertex][0],bounding_box_vertices[vertex][1],bounding_box_vertices[vertex][2]);
		}
		glEnd();
	}
	glEnable(GL_LIGHTING);
}

//------------------------------
//Calculate normal vectors for
//a given mesh and store them in 
//normals
//------------------------------
void process_normals(COORD_ARRAY& coord, VERTEX_INDEX_ARRAY& simplices, int dimension)
{
	if(normals != NULL) //if we already have normals from the last rendered surface delete them
	{
		delete[] normals;
	}
	normals = new Vector[coord.size()/3]; //Total number of vectors = number of vertices = coord.size()/3
	for(int js = 0; js < simplices.size(); js += dimension)
	{
		
		int iv = simplices[js];
		int iv2 = simplices[js + 1];
		int iv3	= simplices[js + 2];	
		Vector v1, v2;
		// n = (p0-p1)X(p1-p2)
		v1[0] = coord[iv*dimension] - coord[iv2*dimension];
		v1[1] = coord[iv*dimension+1] - coord[iv2*dimension+1];
		v1[2] = coord[iv*dimension+2] - coord[iv2*dimension+2];
		v2[0] = coord[iv2*dimension] - coord[iv3*dimension];
		v2[1] = coord[iv2*dimension+1] - coord[iv3*dimension+1];
		v2[2] = coord[iv2*dimension+2] - coord[iv3*dimension+2];
		normals[iv] = normals[iv] + v1.Cross(v2);
		normals[iv2] = normals[iv2] + v1.Cross(v2);
		normals[iv3] = normals[iv3] + v1.Cross(v2);
	}		
	
	for(int i = 0; i < coord.size()/3; i++) //Normalize all vectors
	{
		normals[i].Normalize();
		
	}
	recalculate = false; //We only call this function the first time we render the mesh, subsequent times we want to use the same normal
}

//-------------------------------
//Parse the command line for 
//input file and other parameters
//to be defined later
//-------------------------------
void parse_command_line(int argc, char *argv[])
{
	if(argc > 1)
	{
		io_info.input_filename = argv[1]; //Set input filename
		read_data();
	}
	else
	{
		cerr << "Usage: IsoRender <input_file>\n";
		exit(0);
	}
}

//------------------------------
//Setup io_info to call 
//Marching Cubes Algorithm
//Callback operation of the
//isovalue_spinner GLUI_Spinner
//------------------------------
void setup_marching_cubes(int id)
{
	if(io_info.isovalue.size() > 0) //We only want to generate one isosurface at any given time, make sure there is only one isovalue set
	{
		io_info.isovalue[0] = isovalue; 
	}
	else
	{
		io_info.isovalue.push_back(isovalue);
	}
	io_info.flag_subsample = (subsample == 1) ? false : true; //Set supersample and subsample flags
	io_info.subsample_resolution = subsample;
	io_info.flag_supersample = (supersample == 1) ? false : true;
	io_info.supersample_resolution = supersample;
	if(io_info.flag_subsample && io_info.flag_supersample)
	{
		cout << "Error: Unable to subsample and supersample. Defaulting to subsample." << endl;
		io_info.flag_supersample = false;
	}
	create_isosurface(); //Construct isosurface
}

//---------------------------
//Read in the nrrd file and 
//isotable
//---------------------------
void read_data()
{	
  	IO_TIME io_time = {0.0, 0.0, 0.0};
	IJK::ERROR error;
	MCUBE_TIME mcube_time;
	std::string isotable_directory;
    get_isotable_directory(isotable_directory);
    io_info.isotable_directory = isotable_directory;
	read_nrrd_file(io_info.input_filename, full_scalar_grid,  nrrd_info, io_time);
	if (!check_input(io_info, full_scalar_grid, error)) 
    { 
			throw(error); 
	}
	// Note: mc_data.SetScalarGrid must be called before set_mc_data.
    mc_data.SetScalarGrid(full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
       		io_info.flag_supersample, io_info.supersample_resolution);
    set_mc_data(io_info, mc_data, mcube_time);

    // Note: All flags in mc_data should be set before calling 
    //       read isosurface lookup tables.
    // read isosurface lookup tables
    read_poly_isotable(isotable_directory, mc_data, io_time);
}

//-----------------------------
//Setup parameters to call
//construct_isosurface
//-----------------------------
void create_isosurface()
{
	MCUBE_TIME mcube_time;
  	IO_TIME io_time = {0.0, 0.0, 0.0};
  	IJK::ERROR error;
  	try {

    	std::set_new_handler(memory_exhaustion);

    	// copy nrrd_info into io_info
    	set_io_info(nrrd_info, io_info);


    	// Note: mc_data.SetScalarGrid must be called before set_mc_data.
    	mc_data.SetScalarGrid(full_scalar_grid, io_info.flag_subsample, io_info.subsample_resolution,
       		io_info.flag_supersample, io_info.supersample_resolution);
    	set_mc_data(io_info, mc_data, mcube_time);

    		
		//Set spinner min and max values to that of the scalar grid
		isovalue_spinner->set_float_limits(full_scalar_grid.FindMinScalar(), full_scalar_grid.FindMaxScalar());

    	// construct isosurface or interval volume
    	if (mc_data.IntervalVolumeFlag())
		{
      		construct_interval_volume(io_info, mc_data, mcube_time, io_time);
    	}
    	else 
		{
      		construct_isosurface(io_info, mc_data, mcube_time, io_time);
    	}
	} 
  	catch (ERROR error) 
	{
    	if (error.NumMessages() == 0)
		{	
      		cerr << "Unknown error." << endl;
    	}
    	else { error.Print(cerr); }
    	cerr << "Exiting." << endl;
    	exit(20);
  	}
  	catch (...) 
	{
    	cerr << "Unknown error." << endl;
    	exit(50);
  	};
}

//-----------------------------
//Call Marching Cubes Algorithm
//and create a MC_ISOSURFACE 
//object
//-----------------------------
void construct_isosurface(const IO_INFO & info, const MC_DATA & mc_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
	const int dimension = mc_data.ScalarGrid().Dimension();
  	const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
  	const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

  	io_time.write_time = 0;
  	for (int i = 0; i < info.isovalue.size(); i++) //These loops are not required since there is only one isovalue
	{

    		const SCALAR_TYPE isovalue = info.isovalue[i];

    		MC_ISOSURFACE mc_isosurface;
   		MCUBE_INFO mcube_info(dimension);
    		mcube_info.grid.num_cubes = num_cubes;
    		SNAP_INFO snap_info(dimension);
    		snap_info.grid.num_cubes = num_cubes;

    		if (mc_data.Snap()) 
			{

      			if (info.use_list) 
				{
				std::vector<VERTEX_INDEX> cube_list;
	
				float preprocessing_time;
				IJKSNAPMC::get_nonempty_snap_cubes(mc_data.ScalarGrid(), mc_data.isotable.cube_nep, 
			   		isovalue, cube_list, preprocessing_time);

				mcube_time.preprocessing += preprocessing_time;
				mcube_time.total += preprocessing_time;
	      
				snapMC(mc_data, isovalue, cube_list, mc_isosurface, snap_info);
		      	}
		      	else 
				{
				snapMC(mc_data, isovalue, mc_isosurface, snap_info);
		      	}

		      	mcube_time.Add(snap_info.time);
			}
    		else 
			{
		      marching_cubes(mc_data, isovalue, mc_isosurface, mcube_info);
		      mcube_time.Add(mcube_info.time);
	    	}


   	 	VERTEX_INDEX nums = mc_isosurface.simplex_vert.size()/numv_per_simplex;

    		int grow_factor = 1;
    		int shrink_factor = 1;
    		if (info.flag_subsample) 
      		{ grow_factor = info.subsample_resolution; }
    		if (info.flag_supersample) 
      		{ shrink_factor = info.supersample_resolution; }

    		rescale_vertex_coord(grow_factor, shrink_factor, info.grid_spacing,mc_isosurface.vertex_coord);
		recalculate = true;
		dim = dimension;
    		isosurface.simplex_vert.swap(mc_isosurface.simplex_vert);
    		isosurface.vertex_coord.swap(mc_isosurface.vertex_coord);
		isosurface.cube_containing_simplex.swap(mc_isosurface.cube_containing_simplex);
  	}
}

//----------------------------
//construct_interval_volume
//This operation is never actually used
//----------------------------
void construct_interval_volume(const IO_INFO & info, const MC_DATA & mc_data, MCUBE_TIME & mcube_time, IO_TIME & io_time)
{
	const int dimension = mc_data.ScalarGrid().Dimension();
	const int numv_per_simplex = mc_data.isotable.cube.NumVerticesPerSimplex();
	const int num_cubes = mc_data.ScalarGrid().ComputeNumCubes();

	io_time.write_time = 0;
	for (int i = 0; i+1 < info.isovalue.size(); i++)
	{

		MC_ISOSURFACE mc_ivol;
		MCUBE_INFO mcube_info(dimension);
		mcube_info.grid.num_cubes = num_cubes;
		SNAP_INFO snap_info(dimension);
		snap_info.grid.num_cubes = num_cubes;

		MCVol(mc_data, info.isovalue[i], info.isovalue[i+1], mc_ivol, mcube_info);

		mcube_time.Add(mcube_info.time);


		VERTEX_INDEX nums = mc_ivol.simplex_vert.size()/numv_per_simplex;

		int grow_factor = 1;
		int shrink_factor = 1;

		if (info.flag_subsample) 
		{ grow_factor = info.subsample_resolution; }
		if (info.flag_supersample) 
		{ shrink_factor = info.supersample_resolution; }

		rescale_vertex_coord(grow_factor, shrink_factor, info.grid_spacing,mc_ivol.vertex_coord);
		recalculate = true;
		dim = dimension;
		isosurface.simplex_vert.swap(mc_ivol.simplex_vert);
    		isosurface.vertex_coord.swap(mc_ivol.vertex_coord);
		isosurface.cube_containing_simplex.swap(mc_ivol.cube_containing_simplex);
		
	}
}

//----------------------------
//Set window_x and window_y to
//the center (x,y) position
//----------------------------
void center_on_screen()
{
	window_x = (glutGet (GLUT_SCREEN_WIDTH) - window_width)/2;
	window_y = (glutGet (GLUT_SCREEN_HEIGHT) - window_height)/2;
}

void set_ambient()
{
	float ambient_light[] = {0.4,0.4,0.4,1}; //Seems to work fine
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_light); //Set global ambient term
}

void set_light_positions()
{
	glLightfv(GL_LIGHT0,GL_POSITION,light_position);
}

void write_obj(int id)
{
	if(io_info.isovalue.size() == 0)//If no surface has been rendered there is nothing to write
	{ return;}

	stringstream str;	
	string file;
	string token = io_info.input_filename;
	int start = token.find_last_of("/")+1; //Find start of current directory file
	int end = token.find_last_of(".") - start; //Find position before extension 
	token = token.substr(start, end); //get name of file
	str << token << ".isov=" << int(io_info.isovalue[0]);  //Get full new file name based on sub/super sampling and isovalue
	if(!io_info.flag_subsample && !io_info.flag_supersample)
	{
		str << ".obj";
	}
	else if(io_info.flag_subsample)
	{
		str << ".sub=" << io_info.subsample_resolution << ".obj";
	}
	else
	{
		str << ".super=" << io_info.supersample_resolution << ".obj";
	}
	str >> file;
	ofstream out(file.c_str());
	for(int i = 0; i < isosurface.vertex_coord.size(); i = i + dim) //Output all vertices
	{
		out << "v " <<isosurface.vertex_coord[i] << " " <<isosurface.vertex_coord[i+1] << " "<< isosurface.vertex_coord[i+2] << endl;
	}
	for(int i = 0; i < isosurface.vertex_coord.size()/3; i++) //output all normals
	{
		out << "vn " << normals[i][0] << " " << normals[i][1] << " "<< normals[i][2] << endl;
	}
	for(int js = 0; js < isosurface.simplex_vert.size(); js = js + dim) //output faces
	{
		out << "f ";
		for(int k = 0; k < dim; k++)
		{
			int iv = isosurface.simplex_vert[js + k];
			out << iv +1<< " "; //Note that OBJ file standard starts vertex numbering at 1, not 0. Thus must add one two each index
		}
		out << endl;
	}
	out.close();
	cout << "Wrote table: " << file << endl;
}

void write_off(int id)
{
	if(io_info.isovalue.size() == 0) //If no surface has been rendered there is nothing to write
	{ return;}

	stringstream str;	
	string file;
	string token = io_info.input_filename;
	int start = token.find_last_of("/")+1; //Find start of current directory file
	int end = token.find_last_of(".") - start; //Find position before extension 
	token = token.substr(start, end); //get name of file
	str << token << ".isov=" << int(io_info.isovalue[0]);  //Get full new file name based on sub/super sampling and isovalue
	if(!io_info.flag_subsample && !io_info.flag_supersample) 
	{
		str << ".off";
	}
	else if(io_info.flag_subsample)
	{
		str << ".sub=" << io_info.subsample_resolution << ".off";
	}
	else
	{
		str << ".super=" << io_info.supersample_resolution << ".off";
	}
	str >> file;
	ofstream out(file.c_str());
	out << "OFF" << endl;
	out << isosurface.vertex_coord.size()/3 << " " << isosurface.simplex_vert.size()/3 << " 0"<<endl;  //output number of vertices, faces, and edges
	for(int i = 0; i < isosurface.vertex_coord.size(); i = i + dim)  //Output all vertices
	{
		out << isosurface.vertex_coord[i] << " " <<isosurface.vertex_coord[i+1] << " "<< isosurface.vertex_coord[i+2] << endl;
	}
	for(int js = 0; js < isosurface.simplex_vert.size(); js = js + dim) //output all faces, in this case all faces are triangles
	{
		out << "3 ";
		for(int k = 0; k < dim; k++)
		{
			int iv = isosurface.simplex_vert[js + k];
			out << iv << " ";
		}
		out << endl;
	}
	out.close();
	cout << "Wrote table: " << file << endl;
}

void write_ppm(int id)
{ //Still buggy, may just remove entirely, need to figure out what is wrong
	
	int texture_width, texture_height;
	if(io_info.isovalue.size() == 0)
	{ return;}

	stringstream str;	
	string file;
	string token = io_info.input_filename;
	int start = token.find_last_of("/")+1;
	int end = token.find_last_of(".") - start;
	token = token.substr(start, end);
	str << token << ".isov=" << int(io_info.isovalue[0]);
	if(!io_info.flag_subsample && !io_info.flag_supersample)
	{
		str << ".ppm";
	}
	else if(io_info.flag_subsample)
	{
		str << ".sub=" << io_info.subsample_resolution << ".ppm";
	}
	else
	{
		str << ".super=" << io_info.supersample_resolution << ".ppm";
	}
	str >> file;
	ofstream out(file.c_str(), ios::binary|ios::out);
        int nx, ny, i, j, k;
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
  	GLubyte *image; /* usigned char, 8 bit */
  	nx = viewport[2]-viewport[0];
  	ny =  viewport[3]+viewport[1];
  	image = new GLubyte[nx * ny *3];
  	glPixelStorei (GL_PACK_ALIGNMENT, 1);
  	glReadBuffer(GL_FRONT); 
  	glReadPixels(0, 0, nx, ny, GL_RGB, GL_UNSIGNED_BYTE,image);
  	out<< "P6\n" << nx << " " << viewport[3] << "\n255\n";
	for(int i = ny-1; i >=0; i--)
	{
		if(i >= viewport[1]) //Prevents the black bar at the bottom of the image
		{
			for(int j = 0; j < nx; j++)
			{
				for(int k = 0; k < 3; k++)
					out << (char)image[3*(i*nx+j)+k];
			}
		}
	}
	delete[] image;	
  	out.close();
	cout << "Wrote image: " << file << endl;
}

void reset(int id)
{
	glLoadIdentity(); //Load identity on model view matrix
	for(int i = 0; i < 16; i++) //Set the transformation matrices to identity
	{
		int j = (i == 0 || i == 5 || i == 10 || i == 15);
		matrix1[i] = j;
		matrix2[i] = j;
	}
	glutPostRedisplay();
}

void reshape(int w, int h)
{
	glViewport(0,0,w,h);
	w_top    =  1.0;
	w_bottom = -1.0;
	w_left   = -(double)w/(double)h;
	w_right  = -w_left;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
    	gluPerspective(view_angle,(double)w/(double)h, z_near, z_far); //Set the new aspect ratio

	glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
		case ESC:
		{
			exit(0);
		}
		break;
		default:
		{
			cout << "Unassigned character : " << key << endl;
		}
		break;
	}
}

void mouse(int button, int state, int x, int y)
{
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);
	mouse_x = x;
	mouse_y = y;

	if (state==GLUT_UP)  //Mouse button was pressed
	{
		switch (button)
		{
			case GLUT_LEFT_BUTTON:   
				mouse_left = false; 
			break;
			case GLUT_MIDDLE_BUTTON: 
				mouse_middle = false; 
			break;
			case GLUT_RIGHT_BUTTON:  
				mouse_right = false; 
			break;
		}
	}
	else //mouse buttons are no longer pressed
	{
		switch (button)
		{
			case GLUT_LEFT_BUTTON:	 
				mouse_left = true; 
			break;
			case GLUT_MIDDLE_BUTTON: 
				mouse_middle = true; 
			break;
			case GLUT_RIGHT_BUTTON:	 
				mouse_right = true; 
			break;
		}
	}
	calc_world_coords(drag_x,drag_y,drag_z,x,y,viewport); //Get world coordinates of mouse click
}

void motion(int x, int y)
{
	bool changed = false;
	int dx = x - mouse_x;
	int dy = y - mouse_y;

	int viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport);

	if (dx==0 && dy==0)
	{ return;}

	if (mouse_middle || (mouse_left && mouse_right))
	{
		glLoadIdentity();
		glTranslatef(0,0,dy*0.01);
		glMultMatrixd(matrix1);

		changed = true;
	}
	else if (mouse_left)
	{
		double ax,ay,az;
		double bx,by,bz;
		double angle;

		ax = dy;
		ay = dx;
		az = 0.0;
		angle = sqrt(ax*ax+ay*ay+az*az)/(double)(viewport[2]+1)*180.0;

		bx = matrix2[0]*ax + matrix2[4]*ay + matrix2[8]*az;
		by = matrix2[1]*ax + matrix2[5]*ay + matrix2[9]*az;
		bz = matrix2[2]*ax + matrix2[6]*ay + matrix2[10]*az;
		glRotatef(angle,bx,by,bz);
	
		changed = true;
	}
	else if (mouse_right)
	{
		double px,py,pz;

		calc_world_coords(px,py,pz,x,y,viewport);

		glLoadIdentity();
		glTranslatef(px-drag_x,py-drag_y,pz-drag_z);
		glMultMatrixd(matrix1);

		drag_x = px;
		drag_y = py;
		drag_z = pz;

		changed = true;
	}

	mouse_x = x;
	mouse_y = y;

	if (changed)
	{	
		get_matrix();
		glutPostRedisplay();
	}
}

void calc_world_coords(double &px,double &py,double &pz, int x, int y, int *viewport)
{
	px = (double)(x-viewport[0])/(double)(viewport[2]);
	py = (double)(y-viewport[1])/(double)(viewport[3]);

	px = w_left + (px)*(w_right-w_left);
	py = w_top  + (py)*(w_bottom-w_top);
	pz = z_near;
}

void get_matrix()
{
	glGetDoublev(GL_MODELVIEW_MATRIX,matrix1);
	inverse(matrix1,matrix2);
}


// /********************************************\
// |From Mesa-2.2\src\glu\project.c
// |Compute the inverse of a 4x4 matrix.  
// |Contributed by scotter@lafn.org
// |Thank you.
// \********************************************/
void inverse(const GLdouble *m, GLdouble *out )
{

/* NB. OpenGL Matrices are COLUMN major. */
#define MAT(m,r,c) (m)[(c)*4+(r)]

/* Here's some shorthand converting standard (row,column) to index. */
#define m11 MAT(m,0,0)
#define m12 MAT(m,0,1)
#define m13 MAT(m,0,2)
#define m14 MAT(m,0,3)
#define m21 MAT(m,1,0)
#define m22 MAT(m,1,1)
#define m23 MAT(m,1,2)
#define m24 MAT(m,1,3)
#define m31 MAT(m,2,0)
#define m32 MAT(m,2,1)
#define m33 MAT(m,2,2)
#define m34 MAT(m,2,3)
#define m41 MAT(m,3,0)
#define m42 MAT(m,3,1)
#define m43 MAT(m,3,2)
#define m44 MAT(m,3,3)

   GLdouble det;
   GLdouble d12, d13, d23, d24, d34, d41;
   GLdouble tmp[16]; /* Allow out == in. */

   /* inverse = adjoint / det. (See linear algebra texts.)*/

   /* pre-compute 2x2 dets for last two rows when computing */
   /* cofactors of first two rows. */
   d12 = (m31*m42-m41*m32);
   d13 = (m31*m43-m41*m33);
   d23 = (m32*m43-m42*m33);
   d24 = (m32*m44-m42*m34);
   d34 = (m33*m44-m43*m34);
   d41 = (m34*m41-m44*m31);

   tmp[0] =  (m22 * d34 - m23 * d24 + m24 * d23);
   tmp[1] = -(m21 * d34 + m23 * d41 + m24 * d13);
   tmp[2] =  (m21 * d24 + m22 * d41 + m24 * d12);
   tmp[3] = -(m21 * d23 - m22 * d13 + m23 * d12);

   /* Compute determinant as early as possible using these cofactors. */
   det = m11 * tmp[0] + m12 * tmp[1] + m13 * tmp[2] + m14 * tmp[3];

   /* Run singularity test. */
   if (det == 0.0) {
      /* printf("invert_matrix: Warning: Singular matrix.\n"); */
//	  memcpy(out,_identity,16*sizeof(double));
   }
   else {
      GLdouble invDet = 1.0 / det;
      /* Compute rest of inverse. */
      tmp[0] *= invDet;
      tmp[1] *= invDet;
      tmp[2] *= invDet;
      tmp[3] *= invDet;

      tmp[4] = -(m12 * d34 - m13 * d24 + m14 * d23) * invDet;
      tmp[5] =  (m11 * d34 + m13 * d41 + m14 * d13) * invDet;
      tmp[6] = -(m11 * d24 + m12 * d41 + m14 * d12) * invDet;
      tmp[7] =  (m11 * d23 - m12 * d13 + m13 * d12) * invDet;

      /* Pre-compute 2x2 dets for first two rows when computing */
      /* cofactors of last two rows. */
      d12 = m11*m22-m21*m12;
      d13 = m11*m23-m21*m13;
      d23 = m12*m23-m22*m13;
      d24 = m12*m24-m22*m14;
      d34 = m13*m24-m23*m14;
      d41 = m14*m21-m24*m11;

      tmp[8] =  (m42 * d34 - m43 * d24 + m44 * d23) * invDet;
      tmp[9] = -(m41 * d34 + m43 * d41 + m44 * d13) * invDet;
      tmp[10] =  (m41 * d24 + m42 * d41 + m44 * d12) * invDet;
      tmp[11] = -(m41 * d23 - m42 * d13 + m43 * d12) * invDet;
      tmp[12] = -(m32 * d34 - m33 * d24 + m34 * d23) * invDet;
      tmp[13] =  (m31 * d34 + m33 * d41 + m34 * d13) * invDet;
      tmp[14] = -(m31 * d24 + m32 * d41 + m34 * d12) * invDet;
      tmp[15] =  (m31 * d23 - m32 * d13 + m33 * d12) * invDet;

      memcpy(out, tmp, 16*sizeof(GLdouble));
   }

#undef m11
#undef m12
#undef m13
#undef m14
#undef m21
#undef m22
#undef m23
#undef m24
#undef m31
#undef m32
#undef m33
#undef m34
#undef m41
#undef m42
#undef m43
#undef m44
#undef MAT
}

//-------------------
//Error message for 
//memory exhaustion
//-------------------
void memory_exhaustion()
{
	cerr << "Error: Out of memory.  Terminating program." << endl;
	exit(10);
}
