// csmp_rhino_gmsh.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "CSMP_GMSH_Interface.h"

#include <iostream>
#include <functional>
#include <gmsh.h>
#include <opennurbs_public.h>
#include <sstream>

//#include "bodyType.h"

using namespace std;

namespace csmp {

GMSH_Interface::GMSH_Interface(string input_file_name)
	{
		input_file_name_ = input_file_name;
		MeshSizeMax_ = 1000.;
		MeshSizeFactor_ = 1;
		Algorithm2D_ = 6;
		Algorithm3D_ = 1;
		Tetrahedra_ = 1;
		NbTetrahedra_ = 1;
		NbTriangles_ = 1;
		NbNodes_ = 1;
		Format_ = 16; // 16:VTK
		Smoothing_ = 0; // Number of smoothing steps applied to the final mesh
		MeshSizeFromCurvature_ = 1;
		MinimumElementsPerTwoPi_ = 60;
		ElementOrder_ = 1;

		// initialize gmsh
		printf("\nGmsh: Initializing...");
		gmsh::initialize();

		// initialise OpenNurbs
		ON::Begin();

	}
	
	bool GMSH_Interface::ON_read_3dm_file()
	{
		printf("Reading model from 3dm file...\n");

		return geometry_model_.Read((input_file_name_+".3dm").c_str());
	}
	
	bool  GMSH_Interface::ON_write_3dm_file()
	{
		printf("Writing model to 3dm file...\n");

		return geometry_model_.Write((input_file_name_+"_rewrite.3dm").c_str());

		//printf("**ERROR** Unable to write file\n");
	}

	int GMSH_Interface::ON_add_layer(const wchar_t* layer_name, int red, int green, int blue, int alpha, bool bVisible)
	{
		// alpha: 0 = opaque, 255 = transparent
		// If layer_name is valid, the layer's index (>=0) is returned. Otherwise, ON_UNSET_INT_INDEX is returned

		int layer_index = geometry_model_.AddLayer(layer_name, ON_Color(red, green, blue, alpha));
		if (layer_index != ON_UNSET_INT_INDEX)
		{
			ON_Layer* layer = const_cast<ON_Layer*>(ON_Layer::Cast(geometry_model_.LayerFromIndex(layer_index).ModelComponent()));
			layer->SetVisible(bVisible);
		}
		return layer_index;
	}

	void GMSH_Interface::ON_get_layer_index_to_name_map()
	{
		ONX_ModelComponentIterator it_layer(geometry_model_, ON_ModelComponent::Type::Layer);
		for (const ON_ModelComponent* model_component = it_layer.FirstComponent(); nullptr != model_component; model_component = it_layer.NextComponent())
		{
			const ON_Layer* layer = ON_Layer::Cast(model_component);
			layer_index_name_[layer->Index()].assign(ON_String(layer->Name()).Array());
		}
	}

	// add model component to the OpenNURBS model (we need a way of letting the user know what layer_index values are valid!)
	void GMSH_Interface::ON_add_geometry_component(int layer_index, ON_Geometry* geometry_object)
	{
		// added this function as a convenience for users who want to create their own geometry objects

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}
	
	void GMSH_Interface::ON_add_point(int layer_index, double x, double y, double z)
	{
		// create point
		ON_Point* geometry_object = new ON_Point(ON_3dPoint(x, y, z));

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}

	void GMSH_Interface::ON_add_arccurve(int layer_index, tuple<double, double, double> P_coordinate, tuple<double, double, double> Q_coordinate, tuple<double, double, double> R_coordinate, double angle_in_radians)
	{
		// here is one way of constructing an arc curve using three points (the alternative is from a plane...)
		ON_ArcCurve* geometry_object = new ON_ArcCurve(ON_Arc(ON_Circle(ON_3dPoint(get<0>(P_coordinate), get<1>(P_coordinate), get<2>(P_coordinate)),
			ON_3dPoint(get<0>(Q_coordinate), get<1>(Q_coordinate), get<2>(Q_coordinate)),
			ON_3dPoint(get<0>(R_coordinate), get<1>(R_coordinate), get<2>(R_coordinate))),
			angle_in_radians));

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}

	void GMSH_Interface::ON_add_linecurve(int layer_index, double x_start, double y_start, double z_start, double x_end, double y_end, double z_end)
	{
		// create line curve
		ON_LineCurve* geometry_object = new ON_LineCurve(ON_Line(ON_3dPoint(x_start, y_start, z_start), ON_3dPoint(x_end, y_end, z_end)));

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}

	void GMSH_Interface::ON_add_nurbscurve(int layer_index, int dimension, bool bIsRational, int order, int cv_count, vector<tuple<double, double, double, double>> cv_homogeneous_coordinates, vector<double> knots)
	{
		// check inputs
    // TODO: looks like an error
		if (!knots.size() == order + cv_count - 2 )
		{
			printf("nurbs curve requires the number of knots = order + cv_count - 2\n");
			return;
		}

		// create curve
		ON_NurbsCurve* geometry_object = new ON_NurbsCurve(dimension, bIsRational, order, cv_count);

		for (int idx = 0; idx < geometry_object->CVCount(); idx++)
		{
			double x = get<0>(cv_homogeneous_coordinates[idx]);
			double y = get<1>(cv_homogeneous_coordinates[idx]);
			double z = get<2>(cv_homogeneous_coordinates[idx]);
			double w = get<3>(cv_homogeneous_coordinates[idx]);
			geometry_object->SetCV(idx, ON_4dPoint(x, y, z, w));
		}
		for (int knot_index = 0; knot_index < knots.size(); knot_index++)
			geometry_object->SetKnot(knot_index, knots[knot_index]);

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}

	void GMSH_Interface::ON_add_polycurve(int layer_index, int count)
	{
		printf("add polycurve not implemented as there are too many possibilities (can have any number of any type of curve) - use OpenNURBS to create the object you want\n");
		return;

		// create point
		ON_PolyCurve* geometry_object = new ON_PolyCurve(count);
		for (int idx = 0; idx < count; idx++)
		{
			//geometry_object->Append(ON_Curve)
		}

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}

	void GMSH_Interface::ON_add_polylinecurve(int layer_index, vector<tuple<double, double, double>> point_coordinates)
	{
		// create polylinecurve
		ON_3dPointArray PointArray3d = ON_3dPointArray(point_coordinates.size());
		for (int idx = 0; idx < point_coordinates.size(); idx++)
		{
			PointArray3d[idx] = ON_3dPoint(get<0>(point_coordinates[idx]), get<1>(point_coordinates[idx]), get<2>(point_coordinates[idx]));
		}
		ON_PolylineCurve* geometry_object = new ON_PolylineCurve(PointArray3d);

		// set layer
		ON_3dmObjectAttributes* attributes = new ON_3dmObjectAttributes();
		attributes->m_layer_index = layer_index;

		// add model geometry component
		geometry_model_.AddManagedModelGeometryComponent(geometry_object, attributes);
	}

	void GMSH_Interface::ON_add_extrusion(int layer_index)
	{
		printf("add extrusion not implement\n");
		return;

		ON_Extrusion* geometry_object = new ON_Extrusion();
	}

	void GMSH_Interface::ON_add_revsurface(int layer_index)
	{
		printf("add revsurface not implement\n");
		return;

		ON_RevSurface* geometry_object = new ON_RevSurface();
	}

	void GMSH_Interface::ON_add_sumsurface(int layer_index)
	{
		printf("add sumsurface not implement\n");
		return;

		ON_SumSurface* geometry_object = new ON_SumSurface();
	}

	void GMSH_Interface::ON_add_brep(int layer_index)
	{
		printf("add brep not implement\n");
		return;

		ON_Brep* geometry_object = new ON_Brep();
	}

	void GMSH_Interface::ON_add_nurbssurface(int layer_index, int dimension, bool bIsRational, int order0, int order1, int cv_count0, int cv_count1, vector<tuple<double, double, double, double>> cv_homogeneous_coordinates, vector<double> u_knots, vector<double> v_knots)
	{
		// check inputs
		if ((!u_knots.size() == order0 + cv_count0 - 2) || (!v_knots.size() == order1 + cv_count1 - 2))
		{
			printf("nurbs curve requires the number of knots = order + cv_count - 2\n");
			return;
		}

		// create surface
		ON_NurbsSurface* geometry_object = new ON_NurbsSurface(dimension, bIsRational, order0, order1, cv_count0, cv_count1);

		for (size_t i = 0; i < geometry_object->CVCount(0); i++)
			for (size_t j = 0; j < geometry_object->CVCount(1); j++)
			{
				int idx = i * geometry_object->CVCount(1) + j;
				double x = get<0>(cv_homogeneous_coordinates[idx]);
				double y = get<1>(cv_homogeneous_coordinates[idx]);
				double z = get<2>(cv_homogeneous_coordinates[idx]);
				double w = get<3>(cv_homogeneous_coordinates[idx]);
				geometry_object->SetCV(i, j, ON_4dPoint(x, y, z, w));
			}

		for (int i = 0; i < geometry_object->KnotCount(0); i++)
			geometry_object->SetKnot(0, i, u_knots[i]);

		for (int j = 0; j < geometry_object->KnotCount(1); j++)
			geometry_object->SetKnot(1, j, v_knots[j]);

	}

	void GMSH_Interface::ON_add_planesurface(int layer_index)
	{
		printf("Gmsh does not support infinite planes. Therefore this is not implemented. Use add brep instead to add a trimmed plane surface\n");
		return;

		ON_PlaneSurface geometry_object = ON_PlaneSurface();
	}

	bool GMSH_Interface::AddSurfaceFromCornerPoints(string surface_name, vector<Point<3U>>& points, double mesh_size)
	{
		if (points.size() < 3)
		{
			printf("\n\nAddSurfaceFromCornerPoints: at least 3 points are needed.");
			return false;
		}
		vector<int> points_tags, lines_tags, physical_group_tags;

		// add points
		printf("\nadding points...");
		if (mesh_size == 0)
			mesh_size = points[0].DistanceTo(points[1]) / 2;
		for (unsigned int i = 0; i < points.size(); i++)
		{
			if (mesh_size > points[i].DistanceTo(points[i + 1]) / 2)
				mesh_size = points[i].DistanceTo(points[i + 1]) / 2;

			int point_tag(gmsh::model::occ::addPoint(points[i][0], points[i][1], points[i][2], mesh_size));
			points_tags.push_back(point_tag);
			//physical_group_tags.push_back(point_tag);
			//toolDimTags_.push_back({ 0, point_tag });
		}

		// closing the loop
		points_tags.push_back(points_tags[0]);

		// add lines
		for (unsigned int i = 0; i < points.size(); i++)
		{
			int line_tag(gmsh::model::occ::addLine(points_tags[i], points_tags[i+1]));
			lines_tags.push_back(line_tag);
			//physical_group_tags.push_back(line_tag);
			//toolDimTags_.push_back({ 1, line_tag });
			printf("\nadding lines...%d", line_tag);
		}

		// add curve loop
		printf("\nadding curve loop...");
		int curve_loop_tag(gmsh::model::occ::addCurveLoop(lines_tags));
		printf("%d", curve_loop_tag);

		// add surface
		vector<int>curve_loop_tags;
		curve_loop_tags.push_back(curve_loop_tag);
		printf("\nadding surface...");
		//int surface_tag(gmsh::model::occ::addPlaneSurface(curve_loop_tags));
		int surface_tag(gmsh::model::occ::addSurfaceFilling(curve_loop_tag));
		toolDimTags_.push_back({ 2, surface_tag});
		printf("%d", surface_tag);

		gmsh::model::setEntityName(2, surface_tag, surface_name);

		//physical_group_tags.push_back(surface_tag);

		//surface_tag = gmsh::model::addPhysicalGroup(2, physical_group_tags);

		//gmsh::model::setPhysicalName(2, surface_tag, surface_name);

		bool layer_exist(false);
		int layer_index(0);
		for (int i = 0; i < layers_.size(); i++)
		{
			printf("\nlayer name: %s", layers_[i].name.c_str() );
			if (layers_[i].name == surface_name)
			{
				layer_exist = true;
				layer_index = i;
				break;
			}
		}
		if (layer_exist)
		{
			printf("\nadding surface to an existing layer: %s", layers_[layer_index].name.c_str() );
			layers_[layer_index].objectCount += 1;
		}
		else
		{
			layerStruct new_surface;
			new_surface.name = surface_name;
			new_surface.objectCount = 1;
			new_surface.red = 1;
			new_surface.blue = 10;
			new_surface.green = 20;
			new_surface.alpha = 10;

			printf("\nadding surface to a new layer: %s", new_surface.name.c_str() );
			layers_.push_back(new_surface);
		}
		return false;
	}




	void GMSH_Interface::SetGmshOptions()
	{
		// set options to gaurd against certain things (and take care of possibility we have a function somewhere else that changes them!)

		printf("Gmsh: Setting options that avoid problems...\n");
		gmsh::option::setNumber("General.Verbosity", 5);                    // this is the default verbosity level. It might be useful to increase if we run into problems.
		gmsh::option::setNumber("Geometry.OCCAutoFix", 1);                  // when creating volumes, this ensures that our first surface loop is the exterior, rather than an infinite volume with the surface loop defining a cavity!
		gmsh::option::setNumber("Geometry.OCCBooleanPreserveNumbering", 1); // setting this to zero might be "easier" for Gmsh, but if so then fragmentation will remove entities that aren't modified due to setting removeObject and removeTool, so we set as 1, which is the default anyway

		double tolerance (geometry_model_.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance);
		gmsh::option::setNumber("Geometry.Tolerance", tolerance);
		gmsh::option::setNumber("Geometry.ToleranceBoolean", 2 * tolerance);

		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Setting mesh generation control...\n");
		gmsh::option::setNumber("Mesh.Algorithm", Algorithm2D_);            // Set 2D mesh algorithm: Delaunay (this is the default)
		gmsh::option::setNumber("Mesh.Algorithm3D", Algorithm3D_);          // Set 3D mesh algorithm: Frontal lobe (this is the default)
		gmsh::option::setNumber("Mesh.ElementOrder", ElementOrder_);         // Set element order
		gmsh::option::setNumber("Mesh.MeshSizeMax", MeshSizeMax_);        // Set maximum mesh size (for all the ways to control mesh size see https://gmsh.info//doc/texinfo/gmsh.html#Specifying-mesh-element-sizes)
		gmsh::option::setNumber("Mesh.Format", Format_);					// Mesh output format
		gmsh::option::setNumber("Mesh.MeshSizeFactor", MeshSizeFactor_);
		gmsh::option::setNumber("Mesh.NbNodes", NbNodes_);
		gmsh::option::setNumber("Mesh.NbTetrahedra", NbTetrahedra_);
		gmsh::option::setNumber("Mesh.NbTriangles", NbTriangles_);
		gmsh::option::setNumber("Mesh.Tetrahedra", Tetrahedra_);
		gmsh::option::setNumber("Mesh.Smoothing", Smoothing_);
		gmsh::option::setNumber("Mesh.MeshSizeFromCurvature", MeshSizeFromCurvature_);
		gmsh::option::setNumber("Mesh.MinimumElementsPerTwoPi", MinimumElementsPerTwoPi_);
		gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 2);

	}
	



void GMSH_Interface::Rhino3dm_to_Gmsh(const ON_Geometry* geometry, vector<int>& trimWireTags, gmsh::vectorpair& dimTags)
	{
		// If an application needs to determine if a pointer points to a class derived from ON_SomeClassName, then call ON_SomeClassName::Cast(mystery pointer) and check for a non-null return.
		// ON_Point
		// ON_ArcCurve          not proven
		// ON_CurveOnSurface    I'm not clear what these are. Rhino .NET indicates they are BREP trim curves/loops, and given there is no explicit handling of them in Rhino .NET, which claims to handle all geometry types, I wonder if I need to explicitly handle them here...
		// ON_LineCurve
		// ON_NurbsCurve
		// ON_PolyCurve         not proven
		// ON_PolylineCurve
		// ON_Extrusion         end caps and miter planes not yet implemented
		// ON_NurbsSurface
		// ON_PlaneSurface
		// ON_RevSurface        not proven
		// ON_SumSurface        not proven
		// ON_Brep
		// ...nothing implemented for SubD (Subdivision surface)

		// handle ON_Point
		const ON_Point* point = ON_Point::Cast(geometry);
		if (point != nullptr)
		{
			ON_3dPoint point3d = point->point;
			dimTags.push_back({ 0, gmsh::model::occ::addPoint(point3d.x, point3d.y, point3d.z, 0) });
			return;
		}

		// handle ON_ArcCurve
		const ON_ArcCurve* arc_curve = ON_ArcCurve::Cast(geometry);
		if (arc_curve != nullptr)
		{
			printf("\nON_ArcCurve: Not proven\n");
			ON_Arc arc = arc_curve->m_arc;
			int centerTag = gmsh::model::occ::addPoint(arc.Center().x, arc.Center().y, arc.Center().z, 0);
			if (arc.AngleRadians() < ON_PI)
			{
				int startTag = gmsh::model::occ::addPoint(arc.StartPoint().x, arc.StartPoint().y, arc.StartPoint().z, 0);
				int endTag = gmsh::model::occ::addPoint(arc.EndPoint().x, arc.EndPoint().y, arc.EndPoint().z, 0);
				dimTags.push_back({ 1, gmsh::model::occ::addCircleArc(startTag, centerTag, endTag) });
			}
			else
			{
				ON_Interval interval = arc.DomainRadians();
				double dtheta = arc.AngleRadians() / 3.0;       // since need to be less than pi, this gaurantees it for all all scenarios
				double angle[4] = { interval.m_t[0], interval.m_t[0] + dtheta, interval.m_t[0] + dtheta + dtheta, interval.m_t[1] };
				int endTag;
				for (int idx = 0; idx < 4; idx++)
				{
					int startTag = gmsh::model::occ::addPoint(arc.PointAt(angle[idx]).x, arc.PointAt(angle[idx]).y, arc.PointAt(angle[idx]).z, 0);
					if (idx > 0)
						dimTags.push_back({ 1, gmsh::model::occ::addCircleArc(startTag, centerTag, endTag) });
					endTag = startTag;
				}
			};
			return;
		}

		// handle ON_CurveOnSurface
		const ON_CurveOnSurface* curve_on_surface = ON_CurveOnSurface::Cast(geometry);
		if (curve_on_surface != nullptr)
		{
			printf("\nON_CurveOnSurface: Class handling not implemented!!!\n");
			return;
		}

		// handle ON_LineCurve
		const ON_LineCurve* line_curve = ON_LineCurve::Cast(geometry);
		if (line_curve != nullptr)
		{
			ON_Line line = line_curve->m_line;
			int startTag = gmsh::model::occ::addPoint(line.from.x, line.from.y, line.from.z, 0);
			int endTag = gmsh::model::occ::addPoint(line.to.x, line.to.y, line.to.z, 0);
			dimTags.push_back({ 1, gmsh::model::occ::addLine(startTag, endTag) });
			return;
		}

		// handle ON_NurbsCurve
		const ON_NurbsCurve* nurbs_curve = ON_NurbsCurve::Cast(geometry);
		if (nurbs_curve != nullptr)
		{
			// handle the control vertices
			vector<int> pointTags;
			vector<double> weights;
			for (int cv = 0; cv < nurbs_curve->CVCount(); cv++)
			{
				ON_4dPoint point4d;
				nurbs_curve->GetCV(cv, point4d);
				pointTags.push_back(gmsh::model::occ::addPoint(point4d.x / point4d.w, point4d.y / point4d.w, point4d.z / point4d.w, 0));
				weights.push_back(point4d.w);
			}

			// get the degree
			int degree = nurbs_curve->m_order - 1;

			// handle the knots (openNurbs stores a sorted knot vector with every knot repeated multiplicity times
			// so we need to get the entire list, identify the unique knots, and how many times each is repeated)

			// 1) merge OpenNURBS superfluous knots with OpenNURBS knot vector (both are assumed to be sorted from smallest to largest)
			int mergedknots_n = nurbs_curve->KnotCount() + 2;
			double* mergedknots = (double*)malloc(mergedknots_n * sizeof(double));
			int knot_index1 = 0;
			int knot_index2 = 0;
			for (int knot = 0; knot < mergedknots_n; knot++)
			{
				if (knot_index1 >= nurbs_curve->KnotCount())
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->SuperfluousKnot(knot_index2);
					knot_index2++;
				}
				else if (knot_index2 >= 2)
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->Knot(knot_index1);
					knot_index1++;
				}
				else if (nurbs_curve->Knot(knot_index1) < nurbs_curve->SuperfluousKnot(knot_index2))
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->Knot(knot_index1);
					knot_index1++;
				}
				else
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->SuperfluousKnot(knot_index2);
					knot_index2++;
				}
			}

			// 2) get unique knots and the number of times they occur
			vector<double> knots;
			vector<int> multiplicities;
			for (int knot = 0; knot < mergedknots_n;)
			{
				int knot_sub = knot + 1;
				while (knot_sub < mergedknots_n && mergedknots[knot] == mergedknots[knot_sub])
					knot_sub++;
				knots.push_back(mergedknots[knot]);
				multiplicities.push_back(knot_sub - knot);
				knot = knot_sub;
			}
			free(mergedknots);

			// add B-Spline curve
			dimTags.push_back({ 1, gmsh::model::occ::addBSpline(pointTags, -1, degree, weights, knots, multiplicities) });
			return;
		}

		// handle ON_PolyCurve
		const ON_PolyCurve* poly_curve = ON_PolyCurve::Cast(geometry);
		if (poly_curve != nullptr)
		{
			printf("\nON_PolyCurve: Not proven\n");
			for (int segment_index = 0; poly_curve->Count(); segment_index++)
			{
				Rhino3dm_to_Gmsh(ON_Geometry::Cast(poly_curve->SegmentCurve(segment_index)), trimWireTags, dimTags);
			}
			return;
		}

		// handle ON_PolylineCurve
		const ON_PolylineCurve* polyline_curve = ON_PolylineCurve::Cast(geometry);
		if (polyline_curve != nullptr)
		{
			int endTag;
			for (int point_index = 0; point_index < polyline_curve->PointCount(); point_index++)
			{
				ON_3dPoint point3d = polyline_curve->m_pline[point_index];
				int startTag = gmsh::model::occ::addPoint(point3d.x, point3d.y, point3d.z, 0);
				if (point_index != 0)
				{
					dimTags.push_back({ 1, gmsh::model::occ::addLine(startTag, endTag) });
				}
				endTag = startTag;
			}
			return;
		}

		// handle ON_Extrusion
		const ON_Extrusion* extrusion = ON_Extrusion::Cast(geometry);
		if (extrusion != nullptr)
		{
			printf("\nON_Extrusion: End caps and miter planes at the extremes not yet implemented!!!\n");

			gmsh::vectorpair profileDimTags;
			Rhino3dm_to_Gmsh(ON_Geometry::Cast(extrusion->m_profile), trimWireTags, profileDimTags);

			gmsh::vectorpair outDimTags;
			vector<int> numElements;
			vector<double> heights;
			bool recombine = false;         // this is irrelevant
			gmsh::model::occ::extrude(profileDimTags, extrusion->PathEnd().x - extrusion->PathStart().x, extrusion->PathEnd().y - extrusion->PathStart().y, extrusion->PathEnd().z - extrusion->PathStart().z, outDimTags, numElements, heights, recombine);

			// add end caps using addSurfaceFilling? make volume if both ends capped using add volume obviously
			switch (extrusion->IsCapped())
			{
			case 1: // single bottom cap + sides
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
			case 2: // single top cap + sides
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
			case 3: // bottom and top caps + sides (i.e. a solid)
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
				//vector<int> surfaceTags;
				//for (int idx = 0; idx < outDimTags.size(); idx++)
				//    surfaceTags.push_back(outDimTags[idx].second);
				//vector<int> shellTag;
				//shellTag.push_back(gmsh::model::occ::addSurfaceLoop(surfaceTags, -1, true));
				//outDimTags.clear();
				//outDimTags.push_back({ 3, gmsh::model::occ::addVolume(shellTag,-1) });
				printf("Extrusion is capped, which isn't yet implemented!!!\n");
			}

			// apple miter planes at the extremes using cut somehow?
			switch (extrusion->IsMitered())
			{
			case 1: // start of path is mitered
			case 2: // end of path is mitered
			case 3: // start and end are mitered
				printf("Extrusion is mitered, which isn't yet implemented!!!\n");
			}

			dimTags.insert(dimTags.end(), outDimTags.begin(), outDimTags.end());
			return;
		}

		// handle ON_NurbsSurface
		const ON_NurbsSurface* nurbs_surface = ON_NurbsSurface::Cast(geometry);
		if (nurbs_surface != nullptr)
		{
			// handle the control vertices
			int numPointsU = nurbs_surface->CVCount(0);
			int numPointsV = nurbs_surface->CVCount(1);
			vector<int> pointTags;
			vector<double> weights;
			for (int PointV = 0; PointV < numPointsV; PointV++)
				for (int PointU = 0; PointU < numPointsU; PointU++)
				{
					ON_4dPoint point4d;
					bool error = nurbs_surface->GetCV(PointU, PointV, point4d);
					pointTags.push_back(gmsh::model::occ::addPoint(point4d.x / point4d.w, point4d.y / point4d.w, point4d.z / point4d.w, 0, -1));
					weights.push_back(point4d.w);
				}

			// get the degree
			int degreeU = nurbs_surface->Order(0) - 1;
			int degreeV = nurbs_surface->Order(1) - 1;

			// handle the U knots

			// 1) merge OpenNURBS superfluous knots with OpenNURBS knot vector (both are assumed to be sorted from smallest to largest)
			int mergedknotsU_n = nurbs_surface->KnotCount(0) + 2;
			double* mergedknotsU = (double*)malloc(mergedknotsU_n * sizeof(double));
			int knot_index1 = 0;
			int knot_index2 = 0;
			for (int knot = 0; knot < mergedknotsU_n; knot++)
			{
				if (knot_index1 >= nurbs_surface->KnotCount(0))
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(0, knot_index2);
					knot_index2++;
				}
				else if (knot_index2 >= 2)
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->Knot(0, knot_index1);
					knot_index1++;
				}
				else if (nurbs_surface->Knot(0, knot_index1) < nurbs_surface->SuperfluousKnot(0, knot_index2))
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->Knot(0, knot_index1);
					knot_index1++;
				}
				else
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(0, knot_index2);
					knot_index2++;
				}
			}

			// 2) get unique U knots and the number of times they occur
			vector<double> knotsU;
			vector<int> multiplicitiesU;
			for (int knot = 0; knot < mergedknotsU_n;)
			{
				int knot_sub = knot + 1;
				while (knot_sub < mergedknotsU_n && mergedknotsU[knot] == mergedknotsU[knot_sub])
					knot_sub++;
				knotsU.push_back(mergedknotsU[knot]);
				multiplicitiesU.push_back(knot_sub - knot);
				knot = knot_sub;
			}
			free(mergedknotsU);

			// handle the V knots

			// 1) merge OpenNURBS superfluous knots with OpenNURBS knot vector (both are assumed to be sorted from smallest to largest)
			int mergedknotsV_n = nurbs_surface->KnotCount(1) + 2;
			double* mergedknotsV = (double*)malloc(mergedknotsV_n * sizeof(double));
			knot_index1 = 0;
			knot_index2 = 0;
			for (int knot = 0; knot < mergedknotsV_n; knot++)
			{
				if (knot_index1 >= nurbs_surface->KnotCount(1))
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(1, knot_index2);
					knot_index2++;
				}
				else if (knot_index2 >= 2)
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->Knot(1, knot_index1);
					knot_index1++;
				}
				else if (nurbs_surface->Knot(1, knot_index1) < nurbs_surface->SuperfluousKnot(1, knot_index2))
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->Knot(1, knot_index1);
					knot_index1++;
				}
				else
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(1, knot_index2);
					knot_index2++;
				}
			}

			// 2) get unique V knots and the number of times they occur
			vector<double> knotsV;
			vector<int> multiplicitiesV;
			for (int knot = 0; knot < mergedknotsV_n;)
			{
				int knot_sub = knot + 1;
				while (knot_sub < mergedknotsV_n && mergedknotsV[knot] == mergedknotsV[knot_sub])
					knot_sub++;
				knotsV.push_back(mergedknotsV[knot]);
				multiplicitiesV.push_back(knot_sub - knot);
				knot = knot_sub;
			}
			free(mergedknotsV);

			// add B-Spline surface (sometimes there are OpenCASCADE exceptions if trimWireTags aren't in parametric space)
			dimTags.push_back({ 2, gmsh::model::occ::addBSplineSurface(pointTags, numPointsU, -1, degreeU, degreeV, weights, knotsU, knotsV, multiplicitiesU, multiplicitiesV, trimWireTags, false) });
			return;
		}

		// handle ON_PlaneSurface
		const ON_PlaneSurface* plane_surface = ON_PlaneSurface::Cast(geometry);
		if (plane_surface != nullptr)
		{
			// add plane surface (wire tags must be in 3D space, not parametric space, obviously...)
			dimTags.push_back({ 2, gmsh::model::occ::addPlaneSurface(trimWireTags, -1) });
			return;
		}

		// handle ON_RevSurface
		const ON_RevSurface* rev_surface = ON_RevSurface::Cast(geometry);
		if (rev_surface != nullptr)
		{
			printf("\nON_RevSurface: Not proven\n");
			gmsh::vectorpair reDimTags;
			Rhino3dm_to_Gmsh(ON_Geometry::Cast(rev_surface->m_curve), trimWireTags, reDimTags);

			gmsh::vectorpair outDimTags;
			vector<int> numElements;
			vector<double> heights;
			bool recombine = false;
			gmsh::model::occ::revolve(reDimTags, rev_surface->m_axis.from.x, rev_surface->m_axis.from.y, rev_surface->m_axis.from.z, rev_surface->m_axis.to.x - rev_surface->m_axis.from.x, rev_surface->m_axis.to.y - rev_surface->m_axis.from.y, rev_surface->m_axis.to.z - rev_surface->m_axis.from.z, rev_surface->m_angle.Length(), outDimTags, numElements, heights, recombine);
			dimTags.insert(dimTags.end(), outDimTags.begin(), outDimTags.end());
			return;
		}

		// handle ON_SumSurface
		const ON_SumSurface* sum_surface = ON_SumSurface::Cast(geometry);
		if (sum_surface != nullptr)
		{
			printf("\nON_SumSurface: Not proven\n");

			gmsh::vectorpair curveDimTag;
			Rhino3dm_to_Gmsh(ON_Geometry::Cast(sum_surface->m_curve[0]), trimWireTags, curveDimTag);
			gmsh::vectorpair pathDimTag;
			Rhino3dm_to_Gmsh(ON_Geometry::Cast(sum_surface->m_curve[1]), trimWireTags, pathDimTag);
			vector<int> pathTag;
			pathTag.push_back(pathDimTag[0].second);
			int wireTag = gmsh::model::occ::addWire(pathTag, -1, false);
			gmsh::vectorpair outDimTags;

			// what kind of sweep does Rhino use????  GeomFill_IsCorrectedFrenet is the default for OpenCASCADE, yet Gmsh comments state DiscreteTrihedron seems the most robust; CorrectedFrenet e.g. fails on very simple cases with straight extrusions.
			// I've read on another site that if the extrusion curve has zero curvature points then it will effectively fail, but I presume rhino would use ON_Extrusion for straight extrusion curves, so perhaps we don't have to worry about that?
			// contrary to the Gmsh documentation, the OpenCASCADE will use GeomFill_IsCorrectedFrenet unless the method is GeomFill_IsFrenet, GeomFill_IsCorrectedFrenet, or GeomFill_IsDiscreteTrihedron. I have no idea what these are, but OpenCASCADE documentation
			// states that by default (GeomFill_IsCorrectedFrenet) it Constructs a pipe by sweeping the shape Profile along the wire Spine. The angle made by the spine with the profile is maintained along the length of the pipe. This seems reasonable.
			gmsh::model::occ::addPipe(curveDimTag, wireTag, outDimTags, "CorrectedFrenet");
			dimTags.insert(dimTags.end(), outDimTags.begin(), outDimTags.end());
			return;
		}

		// handle ON_Brep
		const ON_Brep* brep = ON_Brep::Cast(geometry);
		if (brep != nullptr)
		{
			ON_BrepFaceArray brep_face_array = brep->m_F;

			vector<int> trimWireTags;
			gmsh::vectorpair surfaceDimTags;

			for (int face = 0; face < brep_face_array.Count(); face++)
			{
				// we need to know what kind of surface it is, whether trimWireTags is in parametric space or not is dependent on it...
				ON_Surface* surface = brep->m_S[brep_face_array[face].SurfaceIndexOf()];

				trimWireTags.clear();
				if (!brep->FaceIsSurface(face))
					for (int loop = 0; loop < brep_face_array[face].LoopCount(); loop++)
					{
						//ON_BrepLoop* brep_loop = brep_face_array[face].Loop(loop);            // I don't know why this fails, but it does, so I use the next line instead...
						ON_BrepLoop* brep_loop = brep->Loop(brep_face_array[face].m_li[loop]);
						vector<int> curveTags;
						gmsh::vectorpair dimTags;
						for (int trim = 0; trim < brep_loop->TrimCount(); trim++)
						{
							ON_BrepTrim* brep_trim = brep_loop->Trim(trim);

							// not all surfaces support parametric space trim curves in Gmsh, and some surfaces seem to need parametric space trim curves to avoid OpenCASCADE exceptions
							ON_Curve* curve;
							//if (ON_NurbsCurve::Cast(surface) != nullptr)
							if (ON_NurbsSurface::Cast(surface) != nullptr)
							{
								// gmsh nurbs curve can result in OpenCASCADE exceptions if trim wires are in 3D space rather than parametric space
								curve = brep->m_C2[brep_trim->TrimCurveIndexOf()];
							}
							else
							{
								// gmsh plane surface only supports 3D space trim wires
								curve = brep->m_C3[brep_trim->EdgeCurveIndexOf()];
							}
							Rhino3dm_to_Gmsh(ON_Geometry::Cast(curve), trimWireTags, dimTags);   // create a trim curve
							curveTags.push_back(dimTags.back().second);
						}
						trimWireTags.push_back(gmsh::model::occ::addCurveLoop(curveTags, -1));
					}
				Rhino3dm_to_Gmsh(ON_Geometry::Cast(surface), trimWireTags, surfaceDimTags);      // create the trimmed surface
			}

			// make a solid if supposed to be and return a volume tag, otherwise return the surface tags
			if (brep->IsSolid())
			{
				vector<int> surface_Tags;
				for (int face = 0; face < brep_face_array.Count(); face++)
					surface_Tags.push_back(surfaceDimTags[face].second);
				vector<int> shellTags;
				bool sewing = true;
				shellTags.push_back(gmsh::model::occ::addSurfaceLoop(surface_Tags, -1, sewing));
				dimTags.push_back({ 3, gmsh::model::occ::addVolume(shellTags, -1) });
			}
			else
				dimTags.insert(dimTags.end(), surfaceDimTags.begin(), surfaceDimTags.end());

			return;
		}

		// unhandled type!
		printf("Unhandled model geometry class!!!\n");
		return;
	}
 
 
 
 
	void GMSH_Interface::Geometry_to_Gmsh(const ON_Geometry* geometry, vector<int>& trimWireTags, gmsh::vectorpair& dimTags)
	{
		// If an application needs to determine if a pointer points to a class derived from ON_SomeClassName, then call ON_SomeClassName::Cast(mystery pointer) and check for a non-null return.
		// ON_Point
		// ON_ArcCurve
		// ON_CurveOnSurface    I'm not clear what these are. Rhino .NET indicates they are BREP trim curves/loops, and given there is no explicit handling of them in Rhino .NET, which claims to handle all geometry types, I wonder if I need to explicitly handle them here...
		// ON_LineCurve
		// ON_NurbsCurve
		// ON_PolyCurve         not proven
		// ON_PolylineCurve
		// ON_Extrusion         end caps and miter planes not yet implemented
		// ON_NurbsSurface
		// ON_PlaneSurface
		// ON_RevSurface        not proven
		// ON_SumSurface        not proven
		// ON_Brep
		// ...nothing implemented for SubD (Subdivision surface)

		// handle ON_Point
		const ON_Point* point = ON_Point::Cast(geometry);
		if (point != nullptr)
		{
			ON_3dPoint point3d = point->point;
			dimTags.push_back({ 0, gmsh::model::occ::addPoint(point3d.x, point3d.y, point3d.z, 0, -1) });
			return;
		}

		// handle ON_ArcCurve
		const ON_ArcCurve* arc_curve = ON_ArcCurve::Cast(geometry);
		if (arc_curve != nullptr)
		{
			ON_Arc arc = arc_curve->m_arc;
			int centerTag = gmsh::model::occ::addPoint(arc.Center().x, arc.Center().y, arc.Center().z, 0, -1);
			if (arc.AngleRadians() < ON_PI)
			{
				int startTag = gmsh::model::occ::addPoint(arc.StartPoint().x, arc.StartPoint().y, arc.StartPoint().z, 0, -1);
				int endTag = gmsh::model::occ::addPoint(arc.EndPoint().x, arc.EndPoint().y, arc.EndPoint().z, 0, -1);
				dimTags.push_back({ 1, gmsh::model::occ::addCircleArc(startTag, centerTag, endTag, -1) });
			}
			else
			{
				ON_Interval interval = arc.DomainRadians();
				double dtheta = arc.AngleRadians() / 3.0;       // since need to be less than pi, this gaurantees it for all scenarios
				double angle[4] = { interval.m_t[0], interval.m_t[0] + dtheta, interval.m_t[0] + dtheta + dtheta, interval.m_t[1] };
				int endTag;
				for (int idx = 0; idx < 4; idx++)
				{
					int startTag = gmsh::model::occ::addPoint(arc.PointAt(angle[idx]).x, arc.PointAt(angle[idx]).y, arc.PointAt(angle[idx]).z, 0, -1);
					if (idx > 0)
						dimTags.push_back({ 1, gmsh::model::occ::addCircleArc(startTag, centerTag, endTag, -1) });
					endTag = startTag;
				}
			};
			return;
		}

		// handle ON_CurveOnSurface
		const ON_CurveOnSurface* curve_on_surface = ON_CurveOnSurface::Cast(geometry);
		if (curve_on_surface != nullptr)
		{
			printf("\nON_CurveOnSurface: Class handling not implemented!!!\n");
			return;
		}

		// handle ON_LineCurve
		const ON_LineCurve* line_curve = ON_LineCurve::Cast(geometry);
		if (line_curve != nullptr)
		{
			ON_Line line = line_curve->m_line;
			int startTag = gmsh::model::occ::addPoint(line.from.x, line.from.y, line.from.z, 0, -1);
			int endTag = gmsh::model::occ::addPoint(line.to.x, line.to.y, line.to.z, 0, -1);
			dimTags.push_back({ 1, gmsh::model::occ::addLine(startTag, endTag, -1) });
			return;
		}

		// handle ON_NurbsCurve
		const ON_NurbsCurve* nurbs_curve = ON_NurbsCurve::Cast(geometry);
		if (nurbs_curve != nullptr)
		{
			// handle the control vertices
			vector<int> pointTags;
			vector<double> weights;
			for (int cv = 0; cv < nurbs_curve->CVCount(); cv++)
			{
				ON_4dPoint point4d;
				nurbs_curve->GetCV(cv, point4d);
				pointTags.push_back(gmsh::model::occ::addPoint(point4d.x / point4d.w, point4d.y / point4d.w, point4d.z / point4d.w, 0, -1));
				weights.push_back(point4d.w);
			}

			// get the degree
			int degree = nurbs_curve->m_order - 1;

			// handle the knots (openNurbs stores a sorted knot vector with every knot repeated multiplicity times
			// so we need to get the entire list, identify the unique knots, and how many times each is repeated)

			// 1) merge OpenNURBS superfluous knots with OpenNURBS knot vector (both are assumed to be sorted from smallest to largest)
			int mergedknots_n = nurbs_curve->KnotCount() + 2;
			double* mergedknots = (double*)malloc(mergedknots_n * sizeof(double));
			int knot_index1 = 0;
			int knot_index2 = 0;
			for (int knot = 0; knot < mergedknots_n; knot++)
			{
				if (knot_index1 >= nurbs_curve->KnotCount())
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->SuperfluousKnot(knot_index2);
					knot_index2++;
				}
				else if (knot_index2 >= 2)
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->Knot(knot_index1);
					knot_index1++;
				}
				else if (nurbs_curve->Knot(knot_index1) < nurbs_curve->SuperfluousKnot(knot_index2))
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->Knot(knot_index1);
					knot_index1++;
				}
				else
				{
					mergedknots[knot_index1 + knot_index2] = nurbs_curve->SuperfluousKnot(knot_index2);
					knot_index2++;
				}
			}

			// 2) get unique knots and the number of times they occur
			vector<double> knots;
			vector<int> multiplicities;
			for (int knot = 0; knot < mergedknots_n;)
			{
				int knot_sub = knot + 1;
				while (knot_sub < mergedknots_n && mergedknots[knot] == mergedknots[knot_sub])
					knot_sub++;
				knots.push_back(mergedknots[knot]);
				multiplicities.push_back(knot_sub - knot);
				knot = knot_sub;
			}
			free(mergedknots);

			// add B-Spline curve
			dimTags.push_back({ 1, gmsh::model::occ::addBSpline(pointTags, -1, degree, weights, knots, multiplicities) });
			return;
		}

		// handle ON_PolyCurve
		const ON_PolyCurve* poly_curve = ON_PolyCurve::Cast(geometry);
		if (poly_curve != nullptr)
		{
			printf("\nON_PolyCurve: Not proven\n");
			for (int segment_index = 0; poly_curve->Count(); segment_index++)
			{
				Geometry_to_Gmsh(ON_Geometry::Cast(poly_curve->SegmentCurve(segment_index)), trimWireTags, dimTags);
			}
			return;
		}

		// handle ON_PolylineCurve
		const ON_PolylineCurve* polyline_curve = ON_PolylineCurve::Cast(geometry);
		if (polyline_curve != nullptr)
		{
			int endTag;
			for (int point_index = 0; point_index < polyline_curve->PointCount(); point_index++)
			{
				ON_3dPoint point3d = polyline_curve->m_pline[point_index];
				int startTag = gmsh::model::occ::addPoint(point3d.x, point3d.y, point3d.z, 0, -1);
				if (point_index != 0)
				{
					dimTags.push_back({ 1, gmsh::model::occ::addLine(startTag, endTag, -1) });
				}
				endTag = startTag;
			}
			return;
		}

		// handle ON_Extrusion
		const ON_Extrusion* extrusion = ON_Extrusion::Cast(geometry);
		if (extrusion != nullptr)
		{
			printf("\nON_Extrusion: End caps and miter planes at the extremes not yet implemented!!!\n");

			gmsh::vectorpair profileDimTags;
			Geometry_to_Gmsh(ON_Geometry::Cast(extrusion->m_profile), trimWireTags, profileDimTags);

			gmsh::vectorpair outDimTags;
			vector<int> numElements;
			vector<double> heights;
			bool recombine = false;         // this is irrelevant
			gmsh::model::occ::extrude(profileDimTags, extrusion->PathEnd().x - extrusion->PathStart().x, extrusion->PathEnd().y - extrusion->PathStart().y, extrusion->PathEnd().z - extrusion->PathStart().z, outDimTags, numElements, heights, recombine);

			// add end caps using addSurfaceFilling? make volume if both ends capped using add volume obviously
			switch (extrusion->IsCapped())
			{
			case 1: // single bottom cap + sides
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
			case 2: // single top cap + sides
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
			case 3: // bottom and top caps + sides (i.e. a solid)
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
				//outDimTags.push_back({ 2, gmsh::model::occ::addSurfaceFilling(wireTag, -1) });
				//vector<int> surfaceTags;
				//for (int idx = 0; idx < outDimTags.size(); idx++)
				//    surfaceTags.push_back(outDimTags[idx].second);
				//vector<int> shellTag;
				//shellTag.push_back(gmsh::model::occ::addSurfaceLoop(surfaceTags, -1, true));
				//outDimTags.clear();
				//outDimTags.push_back({ 3, gmsh::model::occ::addVolume(shellTag,-1) });
				printf("Extrusion is capped, which isn't yet implemented!!!\n");
			}

			// apple miter planes at the extremes using cut somehow?
			switch (extrusion->IsMitered())
			{
			case 1: // start of path is mitered
			case 2: // end of path is mitered
			case 3: // start and end are mitered
				printf("Extrusion is mitered, which isn't yet implemented!!!\n");
			}

			dimTags.insert(dimTags.end(), outDimTags.begin(), outDimTags.end());
			return;
		}

		// handle ON_NurbsSurface
		const ON_NurbsSurface* nurbs_surface = ON_NurbsSurface::Cast(geometry);
		if (nurbs_surface != nullptr)
		{
			// handle the control vertices
			int numPointsU = nurbs_surface->CVCount(0);
			int numPointsV = nurbs_surface->CVCount(1);
			vector<int> pointTags;
			vector<double> weights;
			for (int PointV = 0; PointV < numPointsV; PointV++)
				for (int PointU = 0; PointU < numPointsU; PointU++)
				{
					ON_4dPoint point4d;
					bool error = nurbs_surface->GetCV(PointU, PointV, point4d);
					pointTags.push_back(gmsh::model::occ::addPoint(point4d.x / point4d.w, point4d.y / point4d.w, point4d.z / point4d.w, 0, -1));
					weights.push_back(point4d.w);
				}

			// get the degree
			int degreeU = nurbs_surface->Order(0) - 1;
			int degreeV = nurbs_surface->Order(1) - 1;

			// handle the U knots

			// 1) merge OpenNURBS superfluous knots with OpenNURBS knot vector (both are assumed to be sorted from smallest to largest)
			int mergedknotsU_n = nurbs_surface->KnotCount(0) + 2;
			double* mergedknotsU = (double*)malloc(mergedknotsU_n * sizeof(double));
			int knot_index1 = 0;
			int knot_index2 = 0;
			for (int knot = 0; knot < mergedknotsU_n; knot++)
			{
				if (knot_index1 >= nurbs_surface->KnotCount(0))
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(0, knot_index2);
					knot_index2++;
				}
				else if (knot_index2 >= 2)
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->Knot(0, knot_index1);
					knot_index1++;
				}
				else if (nurbs_surface->Knot(0, knot_index1) < nurbs_surface->SuperfluousKnot(0, knot_index2))
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->Knot(0, knot_index1);
					knot_index1++;
				}
				else
				{
					mergedknotsU[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(0, knot_index2);
					knot_index2++;
				}
			}

			// 2) get unique U knots and the number of times they occur
			vector<double> knotsU;
			vector<int> multiplicitiesU;
			for (int knot = 0; knot < mergedknotsU_n;)
			{
				int knot_sub = knot + 1;
				while (knot_sub < mergedknotsU_n && mergedknotsU[knot] == mergedknotsU[knot_sub])
					knot_sub++;
				knotsU.push_back(mergedknotsU[knot]);
				multiplicitiesU.push_back(knot_sub - knot);
				knot = knot_sub;
			}
			free(mergedknotsU);

			// handle the V knots

			// 1) merge OpenNURBS superfluous knots with OpenNURBS knot vector (both are assumed to be sorted from smallest to largest)
			int mergedknotsV_n = nurbs_surface->KnotCount(1) + 2;
			double* mergedknotsV = (double*)malloc(mergedknotsV_n * sizeof(double));
			knot_index1 = 0;
			knot_index2 = 0;
			for (int knot = 0; knot < mergedknotsV_n; knot++)
			{
				if (knot_index1 >= nurbs_surface->KnotCount(1))
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(1, knot_index2);
					knot_index2++;
				}
				else if (knot_index2 >= 2)
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->Knot(1, knot_index1);
					knot_index1++;
				}
				else if (nurbs_surface->Knot(1, knot_index1) < nurbs_surface->SuperfluousKnot(1, knot_index2))
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->Knot(1, knot_index1);
					knot_index1++;
				}
				else
				{
					mergedknotsV[knot_index1 + knot_index2] = nurbs_surface->SuperfluousKnot(1, knot_index2);
					knot_index2++;
				}
			}

			// 2) get unique V knots and the number of times they occur
			vector<double> knotsV;
			vector<int> multiplicitiesV;
			for (int knot = 0; knot < mergedknotsV_n;)
			{
				int knot_sub = knot + 1;
				while (knot_sub < mergedknotsV_n && mergedknotsV[knot] == mergedknotsV[knot_sub])
					knot_sub++;
				knotsV.push_back(mergedknotsV[knot]);
				multiplicitiesV.push_back(knot_sub - knot);
				knot = knot_sub;
			}
			free(mergedknotsV);

			// add B-Spline surface (sometimes there are OpenCASCADE exceptions if trimWireTags aren't in parametric space)
			dimTags.push_back({ 2, gmsh::model::occ::addBSplineSurface(pointTags, numPointsU, -1, degreeU, degreeV, weights, knotsU, knotsV, multiplicitiesU, multiplicitiesV, trimWireTags, false) });
			return;
		}

		// handle ON_PlaneSurface
		const ON_PlaneSurface* plane_surface = ON_PlaneSurface::Cast(geometry);
		if (plane_surface != nullptr)
		{
			// add plane surface (wire tags must be in 3D space, not parametric space, obviously...)
			dimTags.push_back({ 2, gmsh::model::occ::addPlaneSurface(trimWireTags, -1) });
			return;
		}

		// handle ON_RevSurface
		const ON_RevSurface* rev_surface = ON_RevSurface::Cast(geometry);
		if (rev_surface != nullptr)
		{
			printf("\nON_RevSurface: Not proven\n");
			gmsh::vectorpair reDimTags;
			Geometry_to_Gmsh(ON_Geometry::Cast(rev_surface->m_curve), trimWireTags, reDimTags);

			gmsh::vectorpair outDimTags;
			vector<int> numElements;
			vector<double> heights;
			bool recombine = false;
			gmsh::model::occ::revolve(reDimTags, rev_surface->m_axis.from.x, rev_surface->m_axis.from.y, rev_surface->m_axis.from.z, rev_surface->m_axis.to.x - rev_surface->m_axis.from.x, rev_surface->m_axis.to.y - rev_surface->m_axis.from.y, rev_surface->m_axis.to.z - rev_surface->m_axis.from.z, rev_surface->m_angle.Length(), outDimTags, numElements, heights, recombine);
			dimTags.insert(dimTags.end(), outDimTags.begin(), outDimTags.end());
			return;
		}

		// handle ON_SumSurface
		const ON_SumSurface* sum_surface = ON_SumSurface::Cast(geometry);
		if (sum_surface != nullptr)
		{
			printf("\nON_SumSurface: Not proven\n");

			gmsh::vectorpair curveDimTag;
			Geometry_to_Gmsh(ON_Geometry::Cast(sum_surface->m_curve[0]), trimWireTags, curveDimTag);
			gmsh::vectorpair pathDimTag;
			Geometry_to_Gmsh(ON_Geometry::Cast(sum_surface->m_curve[1]), trimWireTags, pathDimTag);
			vector<int> pathTag;
			pathTag.push_back(pathDimTag[0].second);
			int wireTag = gmsh::model::occ::addWire(pathTag, -1, false);
			gmsh::vectorpair outDimTags;

			// what kind of sweep does Rhino use????  GeomFill_IsCorrectedFrenet is the default for OpenCASCADE, yet Gmsh comments state DiscreteTrihedron seems the most robust; CorrectedFrenet e.g. fails on very simple cases with straight extrusions.
			// I've read on another site that if the extrusion curve has zero curvature points then it will effectively fail, but I presume rhino would use ON_Extrusion for straight extrusion curves, so perhaps we don't have to worry about that?
			// contrary to the Gmsh documentation, the OpenCASCADE will use GeomFill_IsCorrectedFrenet unless the method is GeomFill_IsFrenet, GeomFill_IsCorrectedFrenet, or GeomFill_IsDiscreteTrihedron. I have no idea what these are, but OpenCASCADE documentation
			// states that by default (GeomFill_IsCorrectedFrenet) it Constructs a pipe by sweeping the shape Profile along the wire Spine. The angle made by the spine with the profile is maintained along the length of the pipe. This seems reasonable.
			gmsh::model::occ::addPipe(curveDimTag, wireTag, outDimTags, "CorrectedFrenet");
			dimTags.insert(dimTags.end(), outDimTags.begin(), outDimTags.end());
			return;
		}

		// handle ON_Brep
		const ON_Brep* brep = ON_Brep::Cast(geometry);
		if (brep != nullptr)
		{
			ON_BrepFaceArray brep_face_array = brep->m_F;

			vector<int> trimWireTags;
			gmsh::vectorpair surfaceDimTags;

			for (int face = 0; face < brep_face_array.Count(); face++)
			{
				// we need to know what kind of surface it is, whether trimWireTags is in parametric space or not is dependent on it...
				ON_Surface* surface = brep->m_S[brep_face_array[face].SurfaceIndexOf()];

				trimWireTags.clear();
				if (!brep->FaceIsSurface(face))
					for (int loop = 0; loop < brep_face_array[face].LoopCount(); loop++)
					{
						//ON_BrepLoop* brep_loop = brep_face_array[face].Loop(loop);            // I don't know why this fails, but it does, so I use the next line instead...
						ON_BrepLoop* brep_loop = brep->Loop(brep_face_array[face].m_li[loop]);
						vector<int> curveTags;
						gmsh::vectorpair dimTags;
						for (int trim = 0; trim < brep_loop->TrimCount(); trim++)
						{
							ON_BrepTrim* brep_trim = brep_loop->Trim(trim);

							// not all surfaces support parametric space trim curves in Gmsh, and some surfaces seem to need parametric space trim curves to avoid OpenCASCADE exceptions
							ON_Curve* curve;
							if (ON_NurbsSurface::Cast(surface) != nullptr)
							{
								// gmsh nurbs curve can result in OpenCASCADE exceptions if trim wires are in 3D space rather than parametric space
								curve = brep->m_C2[brep_trim->TrimCurveIndexOf()];
							}
							else
							{
								// gmsh plane surface only supports 3D space trim wires
								curve = brep->m_C3[brep_trim->EdgeCurveIndexOf()];
							}
							Geometry_to_Gmsh(ON_Geometry::Cast(curve), trimWireTags, dimTags);   // create a trim curve
							curveTags.push_back(dimTags.back().second);
						}
						trimWireTags.push_back(gmsh::model::occ::addCurveLoop(curveTags, -1));
					}
				Geometry_to_Gmsh(ON_Geometry::Cast(surface), trimWireTags, surfaceDimTags);      // create the trimmed surface
			}

			// make a solid if supposed to be and return a volume tag, otherwise return the surface tags
			if (brep->IsSolid())
			{
				vector<int> surfaceTags;
				for (int face = 0; face < brep_face_array.Count(); face++)
					surfaceTags.push_back(surfaceDimTags[face].second);
				vector<int> shellTags;
				bool sewing = true;
				shellTags.push_back(gmsh::model::occ::addSurfaceLoop(surfaceTags, -1, sewing));
				dimTags.push_back({ 3, gmsh::model::occ::addVolume(shellTags, -1) });
			}
			else
				dimTags.insert(dimTags.end(), surfaceDimTags.begin(), surfaceDimTags.end());

			return;
		}

		// unhandled type!
		printf("Unhandled model geometry class!!!\n");
		return;
	}

	// create cracks in the mesh
	void GMSH_Interface::SplitFractures()
	{
		// create cracks in mesh
		gmsh::vectorpair crack_physical_group_dimTags;
		vector<string> split_fractures;
		for (map<string, gmsh::vectorpair>::iterator it = physical_group_name_to_dimTags_.begin(); it != physical_group_name_to_dimTags_.end(); ++it)
		{
//			if (isFractureBody(it->first)) TODO: isFractureBody() is undefined, where is this function?
//			{
				for (int idx = 0; idx < it->second.size(); idx++)
					if (it->second[idx].first == 2)     // just crack surfaces
					{
						crack_physical_group_dimTags.push_back(it->second[idx]);
						split_fractures.push_back(it->first);
					}
//			}
		}
		if (crack_physical_group_dimTags.empty())
			return;
		
		// create cracks in mesh around physical groups as requested by user
		printf("\n--------------------------------------------------------------------------------");
		
		for (int crack_index = 0; crack_index < crack_physical_group_dimTags.size(); crack_index++)
		{
			printf("\nGmsh: Spliting surface: %s\n", split_fractures[crack_index].c_str());
			gmsh::plugin::setNumber("Crack", "Dimension", crack_physical_group_dimTags[crack_index].first);
			gmsh::plugin::setNumber("Crack", "PhysicalGroup", crack_physical_group_dimTags[crack_index].second);
			gmsh::plugin::setNumber("Crack", "OpenBoundaryPhysicalGroup", 0);        // set to >0 if want crack to be left open on that boundary, otherwise lips of the crack sealed
			gmsh::plugin::run("Crack");
		}
	}

	// assess the mesh quality
	void GMSH_Interface::MeshQuality()
	{
		// This is not ideal. It computes statistics for the entire mesh, which include the mesh of objects that we don't transfer to CSMP
		// but were forced to create due to how Gmsh requires us to build objects. Therefore, it is probably better to use CSMP for analysing
		// mesh quality, if it has such a feature available, as we only transfer the relevant parts of the mesh to it.

		// The AnalyseMeshQuality plugin outputs results by printing them to stdout (and the gmsh logger if it is started), but only if General.Verbosity >= 4
		// General.Verbosity (0: silent except for fatal errors, 1: +errors, 2: +warnings, 3: +direct, 4: +information, 5: +status, 99: +debug)
		double old_verbosity;
		gmsh::option::getNumber("General.Verbosity", old_verbosity);
		gmsh::option::setNumber("General.Verbosity", 4);

		// if we want to prevent printing to the screen we could redirect the stdout with something like
		//int stdout_save;
		//char buffer[BIG_ENOUGH];
		//fflush(stdout); //clean everything first
		//stdout_save = dup(STDOUT_FILENO); //save the stdout state
		//freopen("NUL", "a", stdout); //redirect stdout to null pointer
		//setvbuf(stdout, buffer, _IOFBF, 1024); //set buffer to stdout
		//// run plugin
		//freopen("NUL", "a", stdout); //redirect stdout to null again
		//dup2(stdout_save, STDOUT_FILENO); //restore the previous state of stdout
		//setvbuf(stdout, NULL, _IONBF, MDT_BUFSIZE); //disable buffer to print to screen instantly

		// or could do this...
		//#include <sstream>  // obviously whack this up the top!
		//stringstream buffer;
		//streambuf* old = cout.rdbuf(buffer.rdbuf());
		//// run plugin
		//string text = buffer.str();
		//cout.rdbuf(old);

		// plugins print output, so we need to log it!
		gmsh::logger::start();

		// assess mesh quality (option values set below are the defaults. Prints output.)
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Analysing mesh quality...\n");
		gmsh::plugin::setNumber("AnalyseMeshQuality", "JacobianDeterminant", 1);        // compute Jacobian? According to Stephan Matthai this is the only useful thing.
		gmsh::plugin::setNumber("AnalyseMeshQuality", "IGEMeasure", 0);                 // compute IGE?
		gmsh::plugin::setNumber("AnalyseMeshQuality", "ICNMeasure", 0);                 // compute ICN?
		gmsh::plugin::setNumber("AnalyseMeshQuality", "HidingThreshold", 99);
		gmsh::plugin::setNumber("AnalyseMeshQuality", "ThresholdGreater", 1);
		gmsh::plugin::setNumber("AnalyseMeshQuality", "CreateView", 0);                 // create a model-based view of min(J)/max(J), min(IGE) and/or min(ICN)?
		gmsh::plugin::setNumber("AnalyseMeshQuality", "Recompute", 0);                  // force recomputation (set to 1 if the mesh has changed)
		gmsh::plugin::setNumber("AnalyseMeshQuality", "DimensionOfElements", -1);       // analyse elements of the given dimension if equal to 1, 2 or 3; analyse 2D and 3D elements if equal to 4; or analyse elements of the highest dimension if equal to -1.
		gmsh::plugin::run("AnalyseMeshQuality");

		// parse log of plugin printed output
		vector<string> log;
		gmsh::logger::get(log);
		gmsh::logger::stop();
		for (int idx = 0; idx < log.size(); idx++)
		{
			size_t offset = log[idx].find("minJ");
			if (offset != string::npos)
			{
				// what do we do with the results?
				//printf("%s\n", log[idx].data());
			}
		}

		// restore verbosity to original value
		gmsh::option::setNumber("General.Verbosity", old_verbosity);
	}

	// visualise the model and mesh in the Gmsh GUI
	void GMSH_Interface::GmshGUI(bool gui)
	{
		if (gui)
		{
			printf("--------------------------------------------------------------------------------\n");
			printf("Gmsh: Running graphical user interface (close it to contine)...\n");
			gmsh::fltk::run();
		}
	}

	// write CSMP asc and dat files for transferring the Gmsh mesh to CSMP using ICEM element-type names
	void GMSH_Interface::Gmsh_to_ICEM_elementType_string(int gmsh_elementType, string& string_to_append_csmp_element_type)
	{
		// return ICEM element-type names
		switch (gmsh_elementType)                  // write element type
		{                                                                               // GMSH         CSMP                    ICEM
		case 1: string_to_append_csmp_element_type += "BAR_2"; break;                   // MSH_LIN_2    LINEAR_BAR              BAR_2
		case 2: string_to_append_csmp_element_type += "TRI_3"; break;                   // MSH_TRI_3    LINEAR_TRIANGLE         TRI_3
		case 4: string_to_append_csmp_element_type += "TETRA_4"; break;                 // MSH_TET_4    LINEAR_TETRAHEDRON      TETRA_4
		case 8: string_to_append_csmp_element_type += "BAR_3"; break;                   // MSH_LIN_3    QUADRATIC_BAR           BAR_3
		case 9: string_to_append_csmp_element_type += "TRI_6"; break;                   // MSH_TRI_6    QUADRATIC_TRIANGLE      TRI_6
		case 11: string_to_append_csmp_element_type += "TETRA_10"; break;               // MSH_TET_10   QUADRATIC_TETRAHEDRON   TETRA_10
		case 15: string_to_append_csmp_element_type += "UNKNOWN"; break;                // MSH_PNT      POINT_ELEMENT           < no corresponding type >         csp doesn't have this type, which is why we have it as UNKNOWN
		default: string_to_append_csmp_element_type += "UNKNOWN"; break;                // unhandled element type encountered
		}
	}
	
	// write CSMP asc and dat files for transferring the Gmsh mesh to CSMP
	void GMSH_Interface::Gmsh_to_CSMP_elementType_string(int gmsh_elementType, string& string_to_append_csmp_element_type)
	{
		switch (gmsh_elementType)                  // write element type
		{
		case 1: string_to_append_csmp_element_type += "LINEAR_BAR"; break;              // MSH_LIN_2    LINEAR_BAR
		case 2: string_to_append_csmp_element_type += "LINEAR_TRIANGLE"; break;         // MSH_TRI_3    LINEAR_TRIANGLE
		case 4: string_to_append_csmp_element_type += "LINEAR_TETRAHEDRON"; break;      // MSH_TET_4    LINEAR_TETRAHEDRON
		case 8: string_to_append_csmp_element_type += "QUADRATIC_BAR"; break;           // MSH_LIN_3    QUADRATIC_BAR
		case 9: string_to_append_csmp_element_type += "QUADRATIC_TRIANGLE"; break;      // MSH_TRI_6    QUADRATIC_TRIANGLE
		case 11: string_to_append_csmp_element_type += "QUADRATIC_TETRAHEDRON"; break;  // MSH_TET_10   QUADRATIC_TETRAHEDRON
		case 15: string_to_append_csmp_element_type += "UNKNOWN"; break;                // MSH_PNT      POINT_ELEMENT
		default: string_to_append_csmp_element_type += "UNKNOWN"; break;                // unhandled element type encountered
		}
	}

	// return ICEM element-type identifiers
	CSMP_FEM_TYPE GMSH_Interface::Gmsh_to_ICEM_elementType_enum(int gmsh_elementType)
	{
		switch (gmsh_elementType)
		{                                       // GMSH         CSMP                    ICEM
		case 1: return (CSMP_FEM_TYPE)2;         // MSH_LIN_2    LINEAR_BAR              BAR_2 = 2
		case 2: return (CSMP_FEM_TYPE)8;         // MSH_TRI_3    LINEAR_TRIANGLE         TRI_3 = 8
		case 4: return (CSMP_FEM_TYPE)4;         // MSH_TET_4    LINEAR_TETRAHEDRON      TETRA_4 = 4
		case 8: return (CSMP_FEM_TYPE)3;         // MSH_LIN_3    QUADRATIC_BAR           BAR_3 = 3
		case 9: return (CSMP_FEM_TYPE)10;        // MSH_TRI_6    QUADRATIC_TRIANGLE      TRI_6 = 10
		case 11: return (CSMP_FEM_TYPE)5;        // MSH_TET_10   QUADRATIC_TETRAHEDRON   TETRA_10 = 5
		case 15: return UNKNOWN;                // MSH_PNT      POINT_ELEMENT           < no corresponding type >         csp doesn't have this type, which is why we have it as UNKNOWN
		default: return UNKNOWN;                // unhandled element type encountered
		}
	}

	// return CSMP element-type identifiers
	CSMP_FEM_TYPE GMSH_Interface::Gmsh_to_CSMP_elementType_enum(int gmsh_elementType)
	{
		switch (gmsh_elementType)
		{										// GMSH			CSMP  
		case 1: return LINEAR_BAR;              // MSH_LIN_2    LINEAR_BAR
		case 2: return LINEAR_TRIANGLE;         // MSH_TRI_3    LINEAR_TRIANGLE
		case 4: return LINEAR_TETRAHEDRON;      // MSH_TET_4    LINEAR_TETRAHEDRON
		case 8: return QUADRATIC_BAR;           // MSH_LIN_3    QUADRATIC_BAR
		case 9: return QUADRATIC_TRIANGLE;      // MSH_TRI_6    QUADRATIC_TRIANGLE
		case 11: return QUADRATIC_TETRAHEDRON;  // MSH_TET_10   QUADRATIC_TETRAHEDRON
		case 15: return UNKNOWN;                // MSH_PNT      POINT_ELEMENT
		default: return UNKNOWN;                // unhandled element type encountered
		}
	}

	BOX_BOUNDARY GMSH_Interface::Boundary_flag_to_csmp(int flag)
	{
		int outer_flag = flag & 0x3F;
		if (outer_flag != 0)            // use this to distinguish between the bounding box, and something internal to it
			flag = outer_flag;
		switch (flag)
		{
		case 0x00: return BOX_BOUNDARY::NOT;
		case 0x01: return BOX_BOUNDARY::TOP;         // CSP_BOUNDARY_FLAGS::TOP_OUTSIDE;
		case 0x02: return BOX_BOUNDARY::BOTTOM;      // CSP_BOUNDARY_FLAGS::BOTTOM_OUTSIDE;
		case 0x04: return BOX_BOUNDARY::LEFT;        // CSP_BOUNDARY_FLAGS::LEFT_OUTSIDE;
		case 0x05: return BOX_BOUNDARY::EDGE8;       // CSP_BOUNDARY_FLAGS::TOP_LEFT;
		case 0x06: return BOX_BOUNDARY::EDGE5;       // CSP_BOUNDARY_FLAGS::BOTTOM_LEFT;
		case 0x08: return BOX_BOUNDARY::RIGHT;       // CSP_BOUNDARY_FLAGS::RIGHT_OUTSIDE;
		case 0x09: return BOX_BOUNDARY::EDGE7;       // CSP_BOUNDARY_FLAGS::TOP_RIGHT;
		case 0x0A: return BOX_BOUNDARY::EDGE6;       // CSP_BOUNDARY_FLAGS::BOTTOM_RIGHT;
		case 0x10: return BOX_BOUNDARY::FRONT;       // CSP_BOUNDARY_FLAGS::FRONT_OUTSIDE;
		case 0x11: return BOX_BOUNDARY::EDGE11;      // CSP_BOUNDARY_FLAGS::FRONT_TOP;
		case 0x12: return BOX_BOUNDARY::EDGE9;       // CSP_BOUNDARY_FLAGS::FRONT_BOTTOM;
		case 0x14: return BOX_BOUNDARY::EDGE12;      // CSP_BOUNDARY_FLAGS::FRONT_LEFT;
		case 0x15: return BOX_BOUNDARY::CNR8;        // CSP_BOUNDARY_FLAGS::CNR_MAX_MAXZ;          // front top left
		case 0x16: return BOX_BOUNDARY::CNR5;        // CSP_BOUNDARY_FLAGS::CNR_MIN_MAXZ;          // front bottom left
		case 0x18: return BOX_BOUNDARY::EDGE10;      // CSP_BOUNDARY_FLAGS::FRONT_RIGHT;
		case 0x19: return BOX_BOUNDARY::CNR7;        // CSP_BOUNDARY_FLAGS::CNR_MAX;               // front top right
		case 0x1A: return BOX_BOUNDARY::CNR6;        // CSP_BOUNDARY_FLAGS::CNR_MIN_MAXXZ;         // front bottom right
		case 0x20: return BOX_BOUNDARY::BACK;        // CSP_BOUNDARY_FLAGS::BACK_OUTSIDE;
		case 0x21: return BOX_BOUNDARY::EDGE3;       // CSP_BOUNDARY_FLAGS::BACK_TOP;
		case 0x22: return BOX_BOUNDARY::EDGE1;       // CSP_BOUNDARY_FLAGS::BACK_BOTTOM;
		case 0x24: return BOX_BOUNDARY::EDGE4;       // CSP_BOUNDARY_FLAGS::BACK_LEFT;
		case 0x25: return BOX_BOUNDARY::CNR4;        // CSP_BOUNDARY_FLAGS::CNR_MAX_MINXZ;         // back top left
		case 0x26: return BOX_BOUNDARY::CNR1;        // CSP_BOUNDARY_FLAGS::CNR_MIN;               // back bottom left
		case 0x28: return BOX_BOUNDARY::EDGE2;       // CSP_BOUNDARY_FLAGS::BACK_RIGHT;
		case 0x29: return BOX_BOUNDARY::CNR3;        // CSP_BOUNDARY_FLAGS::CNR_MAX_MAXX;          // back top right
		case 0x2A: return BOX_BOUNDARY::CNR2;        // CSP_BOUNDARY_FLAGS::CNR_MIN_MAXX;          // back bottom right
		case 0x40: return BOX_BOUNDARY::INTERNAL;    // CSP_BOUNDARY_FLAGS::REGION_BOUNDARY
		default:
			printf("Unknown flag\n");
			return BOX_BOUNDARY::NOT;
		}
	}

	void GMSH_Interface::OutputMesh()
	{

		// if user cracked mesh, then:  The plugin duplicates the nodes and the elements on the crack and stores them in a new discrete curve (`Dimension' = 1) or surface (`Dimension' = 2).
		// This adds extra nodes and elements to physical group that won't be densely numbered, so we renumber all nodes and elements
		gmsh::option::setNumber("Mesh.SaveAll", 0);                         // This is the default value. It ensures that entities belonging to physical groups have their elements and nodes numbered (assigned tags) first so that it ends up with a dense node and element numbering
		gmsh::model::mesh::renumberNodes();
		gmsh::model::mesh::renumberElements();

		// get mesh nodes and elements (for now we are just grabbing the mesh exactly as Gmsh's MSH file would dump it)
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Writing CSMP model topology (ASC) file...\n");

		// compile the model topology information
		gmsh::vectorpair physicalGroupDimTags;
		gmsh::model::getPhysicalGroups(physicalGroupDimTags, -1);
		map<int, map<string, vector<int>>> model_topology;      // model_topology[elementType][physical group name] = {elementTags}
		size_t nodeTag_max = 0;
		size_t elementTag_max = 0;
		for (int idx = 0; idx < physicalGroupDimTags.size(); idx++)
		{
			string name;
			gmsh::model::getPhysicalName(physicalGroupDimTags[idx].first, physicalGroupDimTags[idx].second, name);

			// Gmsh developer recommendation for an equivalent getElementsForPhysicalGroup to getNodesForPhysicalGroup
			// is getEntitiesForPhysicalGroup then getElements...
			vector<int> tags;
			gmsh::model::getEntitiesForPhysicalGroup(physicalGroupDimTags[idx].first, physicalGroupDimTags[idx].second, tags);
			vector<int> elementTypes;
			vector<vector<size_t>> elementTags;
			vector<vector<size_t>> nodeTags;
			for (int entity_index = 0; entity_index < tags.size(); entity_index++)
			{
				gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTags, physicalGroupDimTags[idx].first, tags[entity_index]);     // get all the elements for the volumes
				for (int elementType_index = 0; elementType_index < elementTypes.size(); elementType_index++)
				{
					model_topology[elementTypes[elementType_index]][name].insert(model_topology[elementTypes[elementType_index]][name].end(), elementTags[elementType_index].begin(), elementTags[elementType_index].end());

					size_t maxTag = *max_element(nodeTags[elementType_index].begin(), nodeTags[elementType_index].end());
					nodeTag_max = max(nodeTag_max, maxTag);
					maxTag = *max_element(elementTags[elementType_index].begin(), elementTags[elementType_index].end());
					elementTag_max = max(elementTag_max, maxTag);
				}
			}
		}

		// create the CSMP model topology file and write the comments
    struct tm ltm;
    time_t now = time(nullptr);

#ifdef _WIN32
    localtime_s(&ltm, &now);
#else
    localtime_r(&now, &ltm);
#endif
		ofstream CSMP_asc(input_file_name_ + ".asc");
		size_t offset = input_file_name_.find_last_of("\\");
		if (offset == string::npos)
			offset = -1;
		CSMP_asc << "Input geometrical model: " << input_file_name_ + ".3dm" <<"\n";                                                                                   // model name
		CSMP_asc << "'" + input_file_name_.substr(offset + 1) + ".asc'" << " generated by CSMP Gmsh interface on " << ltm.tm_mday << "/" << 1 + ltm.tm_mon << "/" << 1900 + ltm.tm_year << "\n";   // comment line (header): 1. filename, 2. date of generation

		// write the CSMP model topology file header
		CSMP_asc << physicalGroupDimTags.size() << " # Number of families\n";
		CSMP_asc << "# Objectname        Elementtype        Material-ID   Number of elements\n";
		for (map<int, map<string, vector<int>>>::iterator it_elementType = model_topology.begin(); it_elementType != model_topology.end(); it_elementType++)
		{
			string asc_line;
			for (map<string, vector<int>>::iterator it = it_elementType->second.begin(); it != it_elementType->second.end(); it++)
			{
				asc_line.assign(it->first);                     // write name
				asc_line.push_back(' ');
				while (asc_line.size() < 20)
					asc_line.push_back(' ');
				Gmsh_to_ICEM_elementType_string(it_elementType->first, asc_line);   // write element type
				asc_line.push_back(' ');
				while (asc_line.size() < 39)
					asc_line.push_back(' ');
				asc_line += "0             ";                   // write material ID (we just set as 0 always)
				asc_line += to_string(it->second.size());  // write number of elements
				CSMP_asc << asc_line << endl;
			}
		}

		// write the CSMP model topology file element ID section
		CSMP_asc << "# now the elements which make up each object are listed in sequence\n";
		for (map<int, map<string, vector<int>>>::iterator it_elementType = model_topology.begin(); it_elementType != model_topology.end(); it_elementType++)
		{
			string asc_line;
			for (map<string, vector<int>>::iterator it = it_elementType->second.begin(); it != it_elementType->second.end(); it++)
			{
				asc_line.assign(it->first);                     // write name
				asc_line.push_back(' ');
				Gmsh_to_ICEM_elementType_string(it_elementType->first, asc_line);   // write element type
				asc_line.push_back(' ');
				asc_line += to_string(it->second.size());  // write number of elements
				asc_line.push_back('\n');

				int idx = 0;
				for (; idx < it->second.size(); idx++)
				{
					if (idx > 0)
						if (idx % 10)                      // write up to 10 element IDs per line
							asc_line.push_back(' ');
						else
							asc_line.push_back('\n');
					asc_line += to_string(it->second[idx] - 1);    // - 1 since element tag numbering starts at 1 in Gmsh (i.e. Mesh.FirstElementTag option default value is 1), whereas element IDs start at 0 in CSMP
				}
				asc_line.push_back('\n');
				CSMP_asc << asc_line;
			}
		}
		CSMP_asc.close();

		// create and write CSMP vset file
		printf("Gmsh: Writing CSMP Vset (BIN) file...\n");
		// https://support.ansys.com/staticassets/ANSYS/Initial%20Content%20Entry/General%20Articles%20-%20Products/ICEM%20CFD%20Interfaces/csp.htm
		// SKUA_CSMP_transfer1.pptx
		ofstream CSMP_dat(input_file_name_ + ".dat", ios::binary);

		vector<size_t> nodeTags;
		vector<double> coord;

		printf("  1. NODE COORDINATE DATA...\n");
		vector<double> px, py, pz;
		for (size_t nodeTag = 1; nodeTag <= nodeTag_max; nodeTag++)
		{
			vector<double> parametricCoord;
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
			gmsh::model::mesh::getNode(nodeTag, coord, parametricCoord,entityDim,entityTag);
			px.push_back(coord[0]);
			py.push_back(coord[1]);
			pz.push_back(coord[2]);
		}
		unsigned long p_size = nodeTag_max;
		CSMP_dat.write((const char*)&p_size, sizeof(p_size));                   // write number of nodes
		CSMP_dat.write((const char*)px.data(), px.size() * sizeof(double));     // write x coordinates
		CSMP_dat.write((const char*)py.data(), py.size() * sizeof(double));     // write y coordinates
		CSMP_dat.write((const char*)pz.data(), pz.size() * sizeof(double));     // write z coordinates

		// 2. NODE FLAGS (this would be more efficient and less prone to me making errors if CSP_BOUNDARY_FLAGS were bit flags...)
		printf("  2. NODE FLAGS...\n");
		map<int, int> node_boundary_flags;
		for (int idx = 0; idx < physicalGroupDimTags.size(); idx++)
		{
			string name;
			gmsh::model::getPhysicalName(physicalGroupDimTags[idx].first, physicalGroupDimTags[idx].second, name);
			gmsh::model::mesh::getNodesForPhysicalGroup(physicalGroupDimTags[idx].first, physicalGroupDimTags[idx].second, nodeTags, coord);
			int flag;
			for (int nodeTags_index = 0; nodeTags_index < nodeTags.size(); nodeTags_index++)
			{
				if (name == "TOP")
					flag = 0x01;
				else if (name == "BOTTOM")
					flag = 0x02;
				else if (name == "LEFT")
					flag = 0x04;
				else if (name == "RIGHT")
					flag = 0x08;
				else if (name == "FRONT")
					flag = 0x10;
				else if (name == "BACK")
					flag = 0x20;
				else
					flag = 0x40;    // all other regions
				node_boundary_flags[nodeTags[nodeTags_index]] |= flag;
			}
		}
		vector<BOX_BOUNDARY> pbflags;
		for (map<int, int>::iterator it = node_boundary_flags.begin(); it != node_boundary_flags.end(); ++it)      // maps are sorted by key (in this case, node tag), which makes our life easy to build up pbflags
		{
			// see CSMP_UserGuide.docx
			pbflags.push_back(Boundary_flag_to_csmp(it->second));
		}
		CSMP_dat.write((const char*)pbflags.data(), pbflags.size() * sizeof(pbflags[0]));    // write node flags

		printf("  3. BOUNDARY CONDITIONS APPLIED TO NODES...\n");
		vector<double> pbvals(pbflags.size(), 0);             // pbvals can all be zeros according to Stephan Matthai.
		CSMP_dat.write((const char*)pbvals.data(), pbvals.size() * sizeof(pbvals[0]));

		printf("  4. TYPE OF FINITE ELEMENTS...\n");
		vector<CSMP_FEM_TYPE> pelmt;
		vector<unsigned long> plist;
		for (size_t elementTag = 1; elementTag <= elementTag_max; elementTag++)
		{
			int elementType;
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
			gmsh::model::mesh::getElement(elementTag, elementType, nodeTags,entityDim,entityTag);
			pelmt.push_back(Gmsh_to_ICEM_elementType_enum(elementType));
			switch (elementType)        // this is to handle when node numbering differs between Gmsh and CSMP
			{
			case 11:                    // second order tetrahedron (MSH_TET_10 in gmsh, QUADRATIC_TETRAHEDRON in csmp)
				size_t temp = nodeTags[8];
				nodeTags[8] = nodeTags[9];
				nodeTags[9] = temp;
			}
			transform(nodeTags.begin(), nodeTags.end(), back_inserter(plist), [](size_t x) { return (unsigned long)x - 1; });   // subtract 1 from all nodes to get from Gmsh base 1 node ID to CSMP base 0 node ID
		}
		unsigned long pelmt_size = pelmt.size();
		CSMP_dat.write((const char*)&pelmt_size, sizeof(pelmt_size));
		CSMP_dat.write((const char*)pelmt.data(), pelmt.size() * sizeof(pelmt[0]));

		printf("  5. MESH CONNECTIVITY...\n");
		unsigned long plist_size = plist.size();
		CSMP_dat.write((const char*)&plist_size, sizeof(plist_size));
		CSMP_dat.write((const char*)plist.data(), plist.size() * sizeof(plist[0]));

		// my general code (somewhat based on: https://gitlab.onelab.info/gmsh/gmsh/blob/master/demos/api/neighbors.py)
		//  1) identify all line vertices, triangle edges, and tetrahedron faces
		map<size_t, vector<pair<size_t, int>>> line_vertices;                                    // line_vertices[vertex 1 node tag] = { elementTag, vertexId }
		map<pair<size_t, size_t>, vector<pair<size_t, int>>> triangle_edges;                // triangle_edges[vertex 1 node tag, vertex 2 node tag] = { elementTag, edgeId }
		map<tuple<size_t, size_t, size_t>, vector<pair<size_t, int >>> tetrahedron_faces;   // tetrahedron_faces[vertex 1 node tag, vertex 2 node tag, vertex 3 node tag] = {elementTag, faceId }
		for (size_t elementTag = 1; elementTag <= elementTag_max; elementTag++)
		{
			int elementType;
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
			gmsh::model::mesh::getElement(elementTag, elementType, nodeTags,entityDim,entityTag);
			switch (elementType)
			{
			case 1:     // 1st order line
			case 8:     // 2nd order line
				line_vertices[nodeTags[1]].push_back({ elementTag, 0 });    // vertex[0]     (opposite node 0)
				line_vertices[nodeTags[0]].push_back({ elementTag, 1 });    // vertex[1]     (opposite node 1)
				break;
			case 2:     // 1st order triangle
			case 9:     // 2nd order triangle
			{
				vector<size_t> edge0 = { nodeTags[1], nodeTags[2] }; sort(edge0.begin(), edge0.end());
				vector<size_t> edge1 = { nodeTags[0], nodeTags[2] }; sort(edge1.begin(), edge1.end());
				vector<size_t> edge2 = { nodeTags[0], nodeTags[1] }; sort(edge2.begin(), edge2.end());
				triangle_edges[{edge0[0], edge0[1]}].push_back({ elementTag, 0 });     // edge[0]  (opposite node 0)
				triangle_edges[{edge1[0], edge1[1]}].push_back({ elementTag, 1 });     // edge[1]  (opposite node 1)
				triangle_edges[{edge2[0], edge2[1]}].push_back({ elementTag, 2 });     // edge[2]  (opposite node 2)
			}
			break;
			case 4:     // 1st order tetrahedron
			case 11:    // 2nd order tetrahedron
			{
				vector<size_t> face0 = { nodeTags[1], nodeTags[2], nodeTags[3] }; sort(face0.begin(), face0.end());
				vector<size_t> face1 = { nodeTags[0], nodeTags[2], nodeTags[3] }; sort(face1.begin(), face1.end());
				vector<size_t> face2 = { nodeTags[0], nodeTags[1], nodeTags[3] }; sort(face2.begin(), face2.end());
				vector<size_t> face3 = { nodeTags[0], nodeTags[1], nodeTags[2] }; sort(face3.begin(), face3.end());
				tetrahedron_faces[{face0[0], face0[1], face0[2]}].push_back({ elementTag, 0 });    // face[0] (opposite node 0)
				tetrahedron_faces[{face1[0], face1[1], face1[2]}].push_back({ elementTag, 1 });    // face[1] (opposite node 1)
				tetrahedron_faces[{face2[0], face2[1], face2[2]}].push_back({ elementTag, 2 });    // face[2] (opposite node 2)
				tetrahedron_faces[{face3[0], face3[1], face3[2]}].push_back({ elementTag, 3 });    // face[3] (opposite node 3)
			}
			break;
			};
		}
		//  2) iterate through each element vertex, edge, and face and identify elements that share that vertex, edge, face
		map<pair<size_t, int>, vector<size_t>> pfverts_temporary;            // pfverts_temporary[elementTag, face/edge/vertex number] = { elementTag that shares this face/edgde/vertex }
		for (map<size_t, vector<pair<size_t, int>>>::iterator it = line_vertices.begin(); it != line_vertices.end(); it++)
		{
			for (int i = 0; i < it->second.size(); i++)
			{
				for (int j = 0; j < it->second.size(); j++)
					if (i != j)
						pfverts_temporary[it->second[i]].push_back(it->second[j].first);
				if (pfverts_temporary[it->second[i]].empty())   // check if there were no neighbours
				{
					if (node_boundary_flags[it->first] == 0)    // if there were no neighbours then it must be on a boundary!
					{
						printf("Line vertex has no neighbour and is not on a boundary!!!\n");
					}
					pfverts_temporary[it->second[i]].push_back(pbflags[it->first - 1]);     // if there are no neighbours then a negative long is used to identify the boundary surface ( -1 since node tags are 1 based, whereas pbflags is 0 based)
				}
			}
		}
		for (map<pair<size_t, size_t>, vector<pair<size_t, int>>>::iterator it = triangle_edges.begin(); it != triangle_edges.end(); it++)
		{
			for (int i = 0; i < it->second.size(); i++)
			{
				for (int j = 0; j < it->second.size(); j++)
					if (i != j)
						pfverts_temporary[it->second[i]].push_back(it->second[j].first);
				if (pfverts_temporary[it->second[i]].empty())
				{
					int flag = node_boundary_flags[it->first.first] & node_boundary_flags[it->first.second];
					if (flag == 0)
					{
						printf("Triangle edge has no neighbour and is not on a boundary!!!\n");
					}
					pfverts_temporary[it->second[i]].push_back(Boundary_flag_to_csmp(flag));
				}
			}
		}
		for (map<tuple<size_t, size_t, size_t>, vector<pair<size_t, int >>>::iterator it = tetrahedron_faces.begin(); it != tetrahedron_faces.end(); it++)
		{
			for (int i = 0; i < it->second.size(); i++)
			{
				for (int j = 0; j < it->second.size(); j++)
					if (i != j)
						pfverts_temporary[it->second[i]].push_back(it->second[j].first);
				if (pfverts_temporary[it->second[i]].empty())
				{
					int flag = node_boundary_flags[get<0>(it->first)] & node_boundary_flags[get<1>(it->first)] & node_boundary_flags[get<2>(it->first)];
					if (flag == 0)
					{
						printf("Tetrahedra face has no neighbour and is not on a boundary!!!\n");
					}
					pfverts_temporary[it->second[i]].push_back(Boundary_flag_to_csmp(flag));
				}
			}
		}
		//  3) potentially more than two elements share the same vertex, edge, or face, but each element should only be assigned one neighbour. Therefore we assign the one that is part of the same physical group, and if there is still more than one neighbour, we chose the one that is best aligned
		vector<long> pfverts;      // holds neighbour information
		for (map<pair<size_t, int>, vector<size_t>>::iterator it = pfverts_temporary.begin(); it != pfverts_temporary.end(); it++)
		{
			if (it->second.empty())
				printf("Element has no type of neighbour/boundary assigned. This should not have happened!!!\n");
			else
			{
				int best_index = 0;
				if (it->second.size() > 1)      // this will only potentially happen for elements with dimension less than 3...
				{
					int elementType;
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
					gmsh::model::mesh::getElement(it->first.first, elementType, nodeTags,entityDim,entityTag);

					// firstly we see if only one of the neighbours is part of the same physical group

					// which physical group(s) are the element part of?
					vector<string> element_groups;
					for (map<string, vector<int>>::iterator it_physical_group = model_topology[elementType].begin(); it_physical_group != model_topology[elementType].end(); ++it_physical_group)
					{
						if (find(it_physical_group->second.begin(), it_physical_group->second.end(), it->first.first) != it_physical_group->second.end())
							element_groups.push_back(it_physical_group->first);
					}

					// which physical group(s) are the neighbouring elements part of?
					vector<vector<string>> neighbour_element_groups(it->second.size());
					for (map<string, vector<int>>::iterator it_physical_group = model_topology[elementType].begin(); it_physical_group != model_topology[elementType].end(); ++it_physical_group)
					{
						for (int idx = 0; idx < it->second.size(); idx++)
							if (find(it_physical_group->second.begin(), it_physical_group->second.end(), it->second[idx]) != it_physical_group->second.end())
								neighbour_element_groups[idx].push_back(it_physical_group->first);
					}

					// which neighbouring elements are in the same physical group as the element?
					vector<string> common;
					vector<int> possible_index;
					for (int idx = 0; idx < it->second.size(); idx++)
					{
						common.clear();
						set_intersection(element_groups.begin(), element_groups.end(), neighbour_element_groups[idx].begin(), neighbour_element_groups[idx].end(), back_inserter(common));
						if (!common.empty())
							possible_index.push_back(idx);
					}

					// was there only one neighbour in the same physical gropu as the element? if so then we are done, otherwise we need to try to further discern based on element and neighbour element alignment
					if (possible_index.size() == 1)
						best_index = possible_index[0];
					else
					{
						// secondly, since there was more than one neighbouring element in the same physical group as the element, we now see which is most aligned

						// if there were neighbouring elements in the same physical group as the element, then lets just keep those
						vector<size_t> temp;
						if (!possible_index.empty())
						{
							for (int idx = 0; idx < possible_index.size(); idx++)
								temp.push_back(it->second[possible_index[idx]]);
							it->second = temp;
						}

						// ok, only elements whose dimension is less than three can potentially have multiple neighbours
						// for line elements we use the dot product to test alignment
						// for triangle elements we use the dot product between vectors normal to the faces
						double metric;
						for (int idx = 0; idx < it->second.size(); idx++)
						{
							vector<size_t> neighbour_nodeTags;
							vector<double> neighbour_coord;
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
							gmsh::model::mesh::getElement(it->second[idx], elementType, neighbour_nodeTags,entityDim,entityTag);

							switch (elementType)
							{
							case 1:     // 1st order line
							case 8:     // 2nd order line
							{
								// find out which is the common vertex and collect the unique nodes in an order that makes it {elementNodeStart, elementNodeEnd/NeighbourElementStart, NeighbourElementEnd}
								vector<size_t> unique_nodes;
								if (nodeTags[0] == neighbour_nodeTags[0])
									unique_nodes = { nodeTags[1], nodeTags[0], neighbour_nodeTags[1] };
								else if (nodeTags[0] == neighbour_nodeTags[1])
									unique_nodes = { nodeTags[1], nodeTags[0], neighbour_nodeTags[0] };
								else if (nodeTags[1] == neighbour_nodeTags[0])
									unique_nodes = { nodeTags[0], nodeTags[1], neighbour_nodeTags[1] };
								else
									unique_nodes = { nodeTags[0], nodeTags[1], neighbour_nodeTags[0] };

								// get the node coordinates
                // TODO: SKM add currently missing object tags; fetch correct values
                int entityDim = 0;
                int entityTag = 0;
								vector<vector<double>> unique_coords;
								for (int index = 0; index < unique_nodes.size(); index++)
								{
									vector<double> parametricCoord;
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
									gmsh::model::mesh::getNode(unique_nodes[index], coord, parametricCoord,entityDim,entityTag);
									unique_coords.push_back(coord);
								}

								// define the element and neighbouring element vectors...
								double v[3] = { unique_coords[1][0] - unique_coords[0][0], unique_coords[1][1] - unique_coords[0][1], unique_coords[1][2] - unique_coords[0][2] };
								double neighbour_v[3] = { unique_coords[2][0] - unique_coords[1][0], unique_coords[2][1] - unique_coords[1][1], unique_coords[2][2] - unique_coords[1][2] };

								// compute the dot product...
								double cosineRadians = (v[0] * neighbour_v[0] + v[1] * neighbour_v[1] + v[2] * neighbour_v[2]) / sqrt((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * (neighbour_v[0] * neighbour_v[0] + neighbour_v[1] * neighbour_v[1] + neighbour_v[2] * neighbour_v[2]));

								// is this result the best so far? (the closer cosineRadians is to +1, the better. Since cosineRadians has a range of -1 to +1, the bigger the better)
								if (idx == 0)
									metric = cosineRadians;
								else
									if (metric < cosineRadians)
									{
										metric = cosineRadians;
										best_index = idx;
									}
							}
							break;
							case 2:     // 1st order triangle
							case 9:     // 2nd order triangle
							{
								// the first three nodes define the face
								nodeTags.erase(nodeTags.begin() + 3, nodeTags.end());
								neighbour_nodeTags.erase(neighbour_nodeTags.begin() + 3, neighbour_nodeTags.end());

								// lets find out which nodes define the common edge and collect the unique nodes in an order {common edge nodeTags, remaining element nodeTag, remaining neighbouring element nodeTag}
								vector<size_t> unique_nodes;
								vector<size_t> edge;
								vector<size_t> remaining_node;
								sort(nodeTags.begin(), nodeTags.end());
								sort(neighbour_nodeTags.begin(), neighbour_nodeTags.end());
								set_intersection(nodeTags.begin(), nodeTags.end(), neighbour_nodeTags.begin(), neighbour_nodeTags.end(), back_inserter(edge));
								unique_nodes.insert(unique_nodes.end(), edge.begin(), edge.end());

								set_difference(nodeTags.begin(), nodeTags.end(), edge.begin(), edge.end(), back_inserter(remaining_node));
								unique_nodes.insert(unique_nodes.end(), remaining_node.begin(), remaining_node.end());

								remaining_node.clear();
								set_difference(neighbour_nodeTags.begin(), neighbour_nodeTags.end(), edge.begin(), edge.end(), back_inserter(remaining_node));
								unique_nodes.insert(unique_nodes.end(), remaining_node.begin(), remaining_node.end());

								// get the node coordinates
// TODO: SKM add currently missing object tags; fetch correct values
int entityDim = 0;
int entityTag = 0;
								vector<vector<double>> unique_coords;
								for (size_t index = 0; index < unique_nodes.size(); index++)
								{
									vector<double> parametricCoord;
									gmsh::model::mesh::getNode( unique_nodes[index], coord, parametricCoord, entityDim, entityTag );
									unique_coords.push_back(coord);
								}

								// compute vectors defining two of the edges of the triangular elements. v0 is the common edge. v1 and v2 are from their respective non-common node to a different common node so as to get the vectors normal to the face we later compute to be in the same direction
								double v0[3] = { unique_coords[1][0] - unique_coords[0][0], unique_coords[1][1] - unique_coords[0][1], unique_coords[1][2] - unique_coords[0][2] };
								double v1[3] = { unique_coords[1][0] - unique_coords[2][0], unique_coords[1][1] - unique_coords[2][1], unique_coords[1][2] - unique_coords[2][2] };
								double v2[3] = { unique_coords[2][0] - unique_coords[0][0], unique_coords[2][1] - unique_coords[0][1], unique_coords[2][2] - unique_coords[0][2] };

								// use cross product to get vectors normal to triangular element faces
								double v[3] = { v0[1] * v1[2] - v0[2] * v1[1], v0[2] * v1[0] - v0[0] * v1[2], v0[0] * v1[1] - v0[1] * v1[0] };
								double neighbour_v[3] = { v0[1] * v2[2] - v0[2] * v2[1], v0[2] * v2[0] - v0[0] * v2[2], v0[0] * v2[1] - v0[1] * v2[0] };

								// compute the dot product...
								double cosineRadians = (v[0] * neighbour_v[0] + v[1] * neighbour_v[1] + v[2] * neighbour_v[2]) / sqrt((v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * (neighbour_v[0] * neighbour_v[0] + neighbour_v[1] * neighbour_v[1] + neighbour_v[2] * neighbour_v[2]));

								// is this result the best so far? (the closer cosineRadians is to +1, the better. Since cosineRadians has a range of -1 to +1, the bigger the better)
								if (idx == 0)
									metric = cosineRadians;
								else
									if (metric < cosineRadians)
									{
										metric = cosineRadians;
										best_index = idx;
									}
							}
							break;
							};
						}

					}

				}
				pfverts.push_back(it->second[best_index]);
			}
		}
		unsigned long pfverts_size = pfverts.size();
		CSMP_dat.write((const char*)&pfverts_size, sizeof(pfverts_size));
		CSMP_dat.write((const char*)pfverts.data(), pfverts.size() * sizeof(pfverts[0]));

		printf("  6. MATERIAL IDENTIFIER FOR EACH ELEMENT = ROCK TYPE...\n");
		vector<unsigned int> pmtrl(elementTag_max, 0);              //  set rock type to 0 for all elements.
		unsigned long pmtrl_size = pmtrl.size();
		CSMP_dat.write((const char*)&pmtrl_size, sizeof(pmtrl_size));
		CSMP_dat.write((const char*)pmtrl.data(), pmtrl.size() * sizeof(pmtrl[0]));

		CSMP_dat.close();
	}

	void GMSH_Interface::AddMeshPoints()
	{
		for (unsigned int i = 0; i < mesh_points_.size(); i++)
		{
			if (mesh_sizes_[i] > 0.)
				gmsh::model::occ::addPoint(mesh_points_[i][0], mesh_points_[i][1], mesh_points_[i][2], mesh_sizes_[i], -1);
			//opennurbs_layer_index_to_gmsh_dimTags_.push_back({ 0, gmsh::model::occ::addPoint(point[0], point[1], point[2], mesh_size, -1) });
		}
	}

	void GMSH_Interface::AddMeshPoints( vector<Point<3U>>& points, vector<double>& mesh_sizes )
	{
		if (mesh_points_.size() == 0)
		{
			mesh_points_ = points;
			mesh_sizes_ = mesh_sizes;
		}
		else
		{
			for (int i = 0; i < points.size(); i++)
			{
				mesh_points_.push_back(points[i]);
				mesh_sizes_.push_back(mesh_sizes[i]);
			}
		}
	};

	bool GMSH_Interface::Run()
	{
		GMSH_Interface::ON_read_3dm_file();

		GMSH_Interface::SetGmshOptions();

		GMSH_Interface::ON_get_layer_index_to_name_map();    // get a map between layer indices and name to identify the desired layer for adding a model component

		/*for (map<int, string>::iterator it = layers.begin(); it != layers.end(); ++it)
		{
			printf("layer[%d] = %s\n", it->first, it->second.data());
			if (it->second.compare("TOP") == 0)
			{
				int layer_index = it->first;
				//opennurbs_model_add_point(layer_index, 0.0, 0.0, 0.0);
				//opennurbs_model_add_arccurve()                    // not properly implemented yet
				//opennurbs_model_add_linecurve();
				//opennurbs_model_add_nurbscurve();
				//opennurbs_model_add_polycurve();                  // not properly implemented yet
				//opennurbs_model_add_polylinecurve();
				//opennurbs_model_add_extrusion();                  // not properly implemented yet
				//opennurbs_model_add_revsurface();                 // not properly implemented yet
				//opennurbs_model_add_sumsurface();                 // not properly implemented yet
				//opennurbs_model_add_brep();                       // not properly implemented yet
				//opennurbs_model_add_nurbssurface();
				//opennurbs_model_add_planesurface();               // not properly implemented yet
			}
		}
		*/
		// write OpenNURBS model to 3dm file
		GMSH_Interface::ON_write_3dm_file();

		// Transferring OpenNURBS model to Gmsh
		GMSH_Interface::ImportGeometry();

		GMSH_Interface::AddMeshPoints();

		// generate mesh	
		GMSH_Interface::Mesh();

		// create cracks in mesh
		GMSH_Interface::SplitFractures();

		// assess mesh quality
		GMSH_Interface::MeshQuality();

		// write CSMP files (argument is the filename base to which the file extension, i.e. either .asc or .dat, will be appended)
		GMSH_Interface::OutputMesh();

		return true;
	}
	
	void GMSH_Interface::ImportGeometry()
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Transferring OpenNURBS model to Gmsh...\n");

		// clear gmsh mode
		gmsh::clear();
		
		// iterate through Rhino 3dm model layers and transfer the model geometry components that are assigned to them
		ONX_ModelComponentIterator it_layer(geometry_model_, ON_ModelComponent::Type::Layer);
		for (const ON_ModelComponent* model_component = it_layer.FirstComponent(); nullptr != model_component; model_component = it_layer.NextComponent())
		{
			const ON_Layer* layer = ON_Layer::Cast(model_component);

			if (layer->IsVisible())  // only transfer visible layers
			{
				// find objects in this Rhino layer and transfer to Gmsh
				ONX_ModelComponentIterator it_modelgeometry(geometry_model_, ON_ModelComponent::Type::ModelGeometry);
				for (const ON_ModelComponent* model_component = it_modelgeometry.FirstComponent(); nullptr != model_component; model_component = it_modelgeometry.NextComponent())
				{
					const ON_ModelGeometryComponent* model_geometry = ON_ModelGeometryComponent::Cast(model_component);
					if (model_geometry->Attributes(nullptr)->m_layer_index == layer->Index())
					{
						vector<int> trimWireTags;      // holds trim wire gmsh tags
						Geometry_to_Gmsh(model_geometry->Geometry(nullptr), trimWireTags, opennurbs_layer_index_to_gmsh_dimTags_[layer->Index()]);       // create Gmsh elementary entities and populate layer_struct.Gmsh_dimTags
					}
				}
				//printf("Layer: %s, objects = %d\n", ON_String(layer->Name()).Array(), (int)opennurbs_layer_index_to_gmsh_dimTags[layer->Index()].size());
			}
		}
	}

	void GMSH_Interface::Mesh()
	{

		// create volume based on bounding surfaces
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Creating bounding volume...\n");

		// construct a vector toolDimTags holding all gmsh entity dim tags that directly map to a opennurbs model geometry component
		// construct a vector holding the bounding surface gmsh tags
		gmsh::vectorpair toolDimTags;
		vector<int> boundingSurfaceTags;
		for (map<int, gmsh::vectorpair>::iterator it = opennurbs_layer_index_to_gmsh_dimTags_.begin(); it != opennurbs_layer_index_to_gmsh_dimTags_.end(); ++it)
		{
			// compile all the dimTags
			toolDimTags.insert(toolDimTags.end(), it->second.begin(), it->second.end());

			// compile all the surface objects in the bounding layers, which will be used to create a bounding volume
			const ON_Layer* layer = ON_Layer::Cast(geometry_model_.LayerFromIndex(it->first).ModelComponent());
			ON_wString layer_name = layer->Name();
			if (!layer_name.CompareOrdinal(L"top", true) ||
				!layer_name.CompareOrdinal(L"bottom", true) ||
				!layer_name.CompareOrdinal(L"left", true) ||
				!layer_name.CompareOrdinal(L"right", true) ||
				!layer_name.CompareOrdinal(L"front", true) ||
				!layer_name.CompareOrdinal(L"back", true))
			{
				for (int index = 0; index < it->second.size(); index++)
					if (it->second[index].first == 2)
						boundingSurfaceTags.push_back(it->second[index].second);
			}
		}


		vector<int> shellTags;
		bool sewing = true;             // when creating surfaces we have no control over boundaries, therefore sewing is required (there is even a comment in the Gmsh source code stating this will be relied upon)
		shellTags.push_back(gmsh::model::occ::addSurfaceLoop(boundingSurfaceTags, -1, sewing));
		gmsh::vectorpair objectDimTags;
		objectDimTags.push_back({ 3, gmsh::model::occ::addVolume(shellTags, -1) });

		// fragment to fragment, embed, and remove duplicates
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Fragmenting, embedding, and removing duplicates (this takes some time!)...\n");
		gmsh::vectorpair outDimTags;
		vector<gmsh::vectorpair> outDimTagsMap;
		bool removeObject = true;
		bool removeTool = true;
		gmsh::model::occ::fragment(objectDimTags, toolDimTags, outDimTags, outDimTagsMap, -1, removeObject, removeTool);

		// synchronise Gmsh with OpenCASCADE model
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Synchronizing with OpenCASCADE model...\n");
		gmsh::model::occ::synchronize();

		// make things look like they do in Rhino (visibility, colour, and physical group)
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Setting entity visibiliy, colour, and physical group as per Rhino3D...\n");
		gmsh::vectorpair dimTags;
		gmsh::model::occ::getEntities(dimTags, -1);         // get the dimTag of every entity
		bool recursive = false;                             // recursive is off since dimTags contains every entity anyway
		gmsh::model::setVisibility(dimTags, 0, recursive);  // turn off visibility of everything (recursive is off since dimTags contains every entity anyway)

		vector<vector<int>> layer_tags;
		gmsh::model::setVisibility(outDimTagsMap[0], 1, recursive);
		int tag;
		for (int idx = 0; idx < outDimTagsMap[0].size(); idx++)
		{
			layer_tags.clear();
			layer_tags.resize(1);
			layer_tags[0].push_back(outDimTagsMap[0][idx].second);
			//printf("\noutDimTagsMap[0][%d].first: %d second: %d", idx, outDimTagsMap[0][idx].first, outDimTagsMap[0][idx].second);
			tag = gmsh::model::addPhysicalGroup(3, layer_tags[0], -1);
			stringstream volume_name;
			volume_name << "MATRIX" << idx;
			gmsh::model::setPhysicalName(3, tag, volume_name.str());
		}

		int start_index = 1;    // index 0 is the volume we created, which we skip over
		for (map<int, gmsh::vectorpair>::iterator it = opennurbs_layer_index_to_gmsh_dimTags_.begin(); it != opennurbs_layer_index_to_gmsh_dimTags_.end(); ++it)
		{
			dimTags.clear();
			layer_tags.clear(); 
			layer_tags.resize(4);
			for (int layer_dimTag_index = 0; layer_dimTag_index < it->second.size(); layer_dimTag_index++)
			{
				dimTags.insert(dimTags.end(), outDimTagsMap[start_index + layer_dimTag_index].begin(), outDimTagsMap[start_index + layer_dimTag_index].end());
			}
			start_index += it->second.size();

			bool recursive = false;                             // recursive is off since dimTags contains every entity anyway
			gmsh::model::setVisibility(dimTags, 1, recursive);  // make the highest dimensional entity visibile, which corresponds to the Rhino object

			const ON_Layer* layer = ON_Layer::Cast(geometry_model_.LayerFromIndex(it->first).ModelComponent());
			gmsh::model::setColor(dimTags, layer->Color().Red(), layer->Color().Green(), layer->Color().Blue(), 255 - layer->Color().Alpha(), recursive);    // set the entity colour (open nurbs 0 = opaque and 255 = transparent, which is opposite of standard definition)
			for (int idx = 0; idx < dimTags.size(); idx++)
			{
				layer_tags[dimTags[idx].first].push_back(dimTags[idx].second);
			}

			// only the highest dimensional entities for each layer should be assigned to a group (i.e. only one physical group
			// per name must be created, i.e. the highest dimension, since ICEM will only import one element type per region)
			string layer_name(ON_String(layer->Name()).Array());
			for (int dim = layer_tags.size() - 1; dim > 0; dim--)   // Saeed doesn't want point (0D) elements, so we start at dim = 1 so as not to create 0D physical groups
			{
				if (layer_tags[dim].size() > 0)
				{
					tag = gmsh::model::addPhysicalGroup(dim, layer_tags[dim]);
					gmsh::model::setPhysicalName(dim, tag, layer_name);
					physical_group_name_to_dimTags_[layer_name].push_back({ dim, tag });
					break;
				}
			}
		}

		// generate mesh
		printf("--------------------------------------------------------------------------------\n");
		printf("Gmsh: Generating mesh (this can take some time!)...\n");
		gmsh::model::mesh::generate(3);        // generate mesh up to 3 dimensions

		gmsh::write(input_file_name_+".msh");

	}

	// finalise OpenNURBS and Gmsh
	GMSH_Interface::~GMSH_Interface()
	{
		printf("--------------------------------------------------------------------------------\n");
		printf("Finalising...\n");

		// end Gmsh
		gmsh::finalize();

		// end OpenNURBS
		ON::End();

		printf("\n\nDONE meshing with Gmsh... ");
	}

}

