template<size_t dim>
void PropertyModeling<dim>::BartonBandisApertureModel() const
{
	const double64 pi = 4.0 * atan(1.0);
	const double64 mu = tan(frac_mech_data_.phi_ / 180. * pi);
	double64 element_aperture, element_perm;
	vector<double64> aperture_list;
	vector<double64> element_unit_norml;
	double64 alpha, theta;
	double64 sigma_n, sigma_s, sigma_eff;
	double64 sigma_n_old, sigma_s_old, sigma_eff_old;
	double64 sigma_nc, sigma_p;
	//const double aj = (frac_mech_data_.JRC0_ * (0.04 * frac_mech_data_.UCS_ / frac_mech_data_.JCS0_ - 0.02)) * 1.e-3; // initial mechanical aperture, Bandis et al. 1983 (m) 
	const double aj = minFractureAperture_; // initial mechanical aperture
	double64 delVj, Vm, Kni;
	double64 deln, delp, delm, M, dilm, phim, phimd;
	double64 Ln, JRCn, JCSn, JRCm, i;
	const double64 A(-0.2960), B(-0.0056), C(2.2410), D(-0.2450); // fitting parameters for maximum normal closure, Bandis et al. 1983
	const double64 E(8.57), F(-0.68);						  // fitting parameters for maximum normal closure, Bandis et al. 1983
	const double64 L0(100.e-3); // lab scale sample which have been analzyed by Bandis et al. 1983
	double64 list[] = { 0., 0.3, 0.6, 1., 2., 4., 10., 25., 100. };
	vector<double64> del_2_delp;
	del_2_delp.assign(list, list + sizeof(list) / sizeof(double64));
	vector<double64> pre_peak, post_peak, JRCm_vec, del_vec;
	size_t pre_peak_count, post_peak_count;
	bool shear_failure;
	size_t peak_idx = find(del_2_delp.begin(), del_2_delp.end(), 1.) - del_2_delp.begin();

	string fileName = modelName_ + "_" + "stress_input.txt";
	ofstream file(fileName.c_str());
	file << "phi              = " << frac_mech_data_.phi_ << " [degrees]" << endl;
	file << "JRC0             = " << frac_mech_data_.JRC0_ << " [-]" << endl;
	file << "JCS0             = " << frac_mech_data_.JCS0_ << " [Pa]" << endl;
	file << "UCS              = " << frac_mech_data_.UCS_ << " [Pa]" << endl;
	file << "initial aperture = " << aj << " [m]" << endl;
	file << "maximum stress   = " << sigma1_ << " [Pa]" << endl;
	file << "minimum stress   = " << sigma3_ << " [Pa]" << endl;
	file << "orientation      = " << beta_ << " [degrees]" << endl;
	file << "fluid pressure   = " << pf_ << " [Pa]" << endl << endl;
	/*
	file << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
	file << setw(10) << "fracture" << setw(12) << "JRCn" << setw(14) << "element no." << setw(10) << "sigma_n" << setw(15) << "sigma_s" << setw(14) << "deln" << setw(16) << "JRCm" << setw(15) << "dilm" << setw(17) << "aperture" << setw(12) << "failure" << endl;
	file << "----------------------------------------------------------------------------------------------------------------------------------------" << endl;
	*/

	size_t counter = 0;

	string modelRegionName;
	const size_t fracNameLength = fracName_.size();
	typename map<string, Region<dim> >::const_iterator it_r;
	for (it_r = model_->UniqueRegionsBegin(); it_r != model_->UniqueRegionsEnd(); ++it_r) {
		modelRegionName = (*it_r).first;
		if (modelRegionName.compare(0, fracNameLength, fracName_.c_str()) == 0) {
			Region<dim>& r_ref(model_->Region(modelRegionName.c_str()));
			if (truncation_censoring_flag_) {
				map<string, double64 >::const_iterator it;
				it = updatedFractureLength_.find(modelRegionName);
				Ln = (*it).second;
			}
			else {
				Ln = r_ref.Volume();
			}
			JRCn = frac_mech_data_.JRC0_ * pow(Ln / L0, -0.02 * frac_mech_data_.JRC0_);
			JCSn = frac_mech_data_.JCS0_ * pow(Ln / L0, -0.03 * frac_mech_data_.JRC0_);
			delp = Ln / 500. * pow(JRCn / Ln, 1. / 3.); // peak shear displacement
			vector<Element<dim>*> parents;
			vector<Point<dim>> element_nodes;
			//cout << "region = " << modelRegionName << "   elements = " << r_ref.Elements() << endl;
			for (typename vector<Element<dim>*>::const_iterator it_e = r_ref.ElementsBegin(); it_e != r_ref.ElementsEnd(); ++it_e) {
				counter += 1;
				deln = 0.;
				delm = 0.;
				dilm = 0.;
				JRCm = 0.;
				shear_failure = false;
				(*it_e)->UnitNormal(element_unit_norml);
				double64 len_norm = sqrt(element_unit_norml[0] * element_unit_norml[0] + element_unit_norml[1] * element_unit_norml[1]);
				double64 len_sigma = sqrt(sigma1_ * cos(beta_ * pi / 180) * sigma1_ * cos(beta_ * pi / 180) + sigma1_ * sin(beta_ * pi / 180) * sigma1_ * sin(beta_ * pi / 180));
				alpha = atan2(element_unit_norml[1], element_unit_norml[0]);
				theta = acos((element_unit_norml[0] * sigma1_ * cos(beta_ * pi / 180) + element_unit_norml[1] * sigma1_ * sin(beta_ * pi / 180)) / (len_norm * len_sigma));
				sigma_n_old = (sigma1_ + sigma3_) / 2 + (sigma1_ - sigma3_) / 2 * cos(2 * theta);
				sigma_s_old = (sigma1_ - sigma3_) / 2 * sin(2 * theta);
				sigma_eff_old = sigma_n_old - pf_;
				sigma_n = sigma1_ * cos(theta) * cos(theta) + sigma3_ * sin(theta) * sin(theta);
				(*(*it_e)).Store(model_->Database().StorageKey("normal stress"), makeScalar(PLAIN, sigma_n));
				sigma_s = (sigma1_ - sigma3_) * sin(theta) * cos(theta);
				(*(*it_e)).Store(model_->Database().StorageKey("shear stress"), makeScalar(PLAIN, sigma_s));
				sigma_eff = sigma_n - pf_;
				(*(*it_e)).Store(model_->Database().StorageKey("effective stress"), makeScalar(PLAIN, sigma_eff));
				//cout << counter << "   " << modelRegionName << "   " << (*it_e)->Idx() << "   " << sigma_eff << endl;				
				if (sigma_eff > 0.) {
					// 1) normal displacement model (Baghbanan and Jing 2008)
					//sigma_nc = (0.487 * aj * 1.e6 + 2.51) * 1.e6; // (Pa)
					//delVj = 9. * aj * sigma_eff / (10. * sigma_eff + sigma_nc); // normal clouser
					// 2) normal displacement model (Bandis et al. 1983)
					Kni = (-7.15 + 1.75 * frac_mech_data_.JRC0_ + 0.02 * (frac_mech_data_.JCS0_ * 1e-6 / (aj * 1e3))) * 1.e9; //  initial nornal stiffness (Pa/m)
					//Vm = (A + B * frac_mech_data_.JRC0_ + C * pow(frac_mech_data_.JCS0_ * 1e-6 / (aj * 1e3), D)) * 1e-3; // maximum normal closure, Bandis et al. 1983
					Vm = (E * pow(frac_mech_data_.JCS0_ * 1e-6 / (aj * 1e3), F)) * 1e-3; // maximum normal closure, Bandis et al. 1983 (m)
					delVj = sigma_eff * Vm / (Kni * Vm + sigma_eff);
					deln = Vm - delVj;
					// shear failure check (condition for frictional sliding)
					if (fabs(sigma_s) > fabs(mu * sigma_eff)) {
						shear_failure = true;
						// shear dilation model (Barton 1982, Olson and Barton 2001)
						sigma_p = sigma_eff * tan((JRCn * log10(JCSn / sigma_eff) + frac_mech_data_.phi_) / 180 * pi); // peak shear strength
						JRCm = (atan(abs(sigma_s) / sigma_eff) - frac_mech_data_.phi_ * pi / 180) / log10(JCSn / sigma_eff);
						i = frac_mech_data_.JRC0_ * log10(frac_mech_data_.JCS0_ / sigma_eff); // initial dilation angle (degrees)
						double64 list[] = { -frac_mech_data_.phi_ / i, 0., 0.75, 1., 0.85, 0.7, 0.5, 0.4, 0. };
						vector<double64> JRCm_2_JRCp;
						JRCm_2_JRCp.assign(list, list + sizeof(list) / sizeof(double64));
						if (JRCm > 0. && JRCm <= JRCn) {
							if (abs(sigma_s) < sigma_p) {
								vector<double64> x_pre_peak, y_pre_peak;
								x_pre_peak.assign(JRCm_2_JRCp.begin(), JRCm_2_JRCp.begin() + peak_idx + 1);
								y_pre_peak.assign(del_2_delp.begin(), del_2_delp.begin() + peak_idx + 1);
								pre_peak_count = count_if(x_pre_peak.begin(), x_pre_peak.end(), bind2nd(less<double>(), JRCm / JRCn));
								Point<1U> p1, p2, pt;
								p1 = x_pre_peak[pre_peak_count];
								double64 p1_val = y_pre_peak[pre_peak_count];
								p2 = x_pre_peak[pre_peak_count - 1];
								double64 p2_val = y_pre_peak[pre_peak_count - 1];
								pt = JRCm / JRCn;
								delm = delp * linearInterpolate(make_pair(p1, p1_val), make_pair(p2, p2_val), pt);
							}
							else {
								vector<double64> x_post_peak, y_post_peak;
								x_post_peak.assign(JRCm_2_JRCp.begin() + peak_idx, JRCm_2_JRCp.end());
								y_post_peak.assign(del_2_delp.begin() + peak_idx, del_2_delp.end());
								post_peak_count = count_if(x_post_peak.begin(), x_post_peak.end(), bind2nd(less<double>(), JRCm / JRCn));
								Point<1U> p1, p2, pt;
								p1 = x_post_peak[x_post_peak.size() - post_peak_count];
								double64 p1_val = y_post_peak[y_post_peak.size() - post_peak_count];
								p2 = x_post_peak[x_post_peak.size() - post_peak_count - 1];
								double64 p2_val = y_post_peak[y_post_peak.size() - post_peak_count - 1];
								pt = JRCm / JRCn;
								delm = delp * linearInterpolate(make_pair(p1, p1_val), make_pair(p2, p2_val), pt);
							}
							M = JRCn / (12 * log10(JCSn / sigma_eff)) + 0.70; // damage coefficient
							phim = JRCm * log10(JCSn / sigma_eff) + frac_mech_data_.phi_; // mobilized friction angle
							phimd = 1 / M * JRCm * log10(JCSn / sigma_eff); // mobilized dialtion angle
							dilm = delm * tan(phimd / 180 * pi);
						}
					}
				}
				element_aperture = aj - deln + dilm;
				element_perm = ParallelPlatePermeability(element_aperture);

				if (element_aperture > maxFractureAperture_) {
					element_aperture = maxFractureAperture_;
					element_perm = ParallelPlatePermeability(element_aperture);
				}

				if (element_aperture < minFractureAperture_) {
					element_aperture = minFractureAperture_;
					if (!constMatrixPermFlag_) {
						// parents which have a common fracture element
						parents.clear();
						element_nodes.clear();
						for (vector<Node<dim>*>::const_iterator it_n = (*it_e)->NodesBegin(); it_n != (*it_e)->NodesEnd(); ++it_n) {
							element_nodes.push_back((*it_n)->Coordinate());
							for (size_t parent = 0; parent < (*it_n)->Parents(); ++parent) {
								parents.push_back((*it_n)->Parent(parent));
							}
						}
						sort(parents.begin(), parents.end());
						set<Element<dim>*> unique_parents(parents.begin(), parents.end());
						std::list<Element<dim>*> common_parents;
						set_difference(parents.begin(), parents.end(),
							unique_parents.begin(), unique_parents.end(),
							back_inserter(common_parents));
						std::list<Element<dim>*>::iterator fracture_element_position = common_parents.begin();
						fracture_element_position = find(common_parents.begin(), common_parents.end(), (*it_e));
						common_parents.erase(fracture_element_position);
						double64 parent_perm_avg(0.);
						for (std::list<Element<dim>*>::const_iterator it_ne = common_parents.begin(); it_ne != common_parents.end(); ++it_ne) {
							element_perm += (*it_ne)->Read(permeability_key_);
						}
						element_perm = element_perm / common_parents.size();
					}
					else {
						element_perm = constMatrixPerm_;
					}
				}

				aperture_list.push_back(element_aperture);
				(*(*it_e)).Store(model_->Database().StorageKey("volume modifier"), makeScalar(PLAIN, element_aperture));
				(*(*it_e)).Store(model_->Database().StorageKey("permeability"), makeScalar(PLAIN, element_perm));

				/*
				file << setw(10) << modelRegionName;
				file << setw(10) << setprecision(3) << JRCn;
				file << setw(10) << (*it_e)->Idx();
				file << setw(15) << setprecision(3) << sigma_n;
				file << setw(15) << setprecision(3) << sigma_s;
				file << setw(15) << setprecision(3) << deln;
				file << setw(15) << setprecision(3) << JRCm;
				file << setw(15) << setprecision(3) << dilm;
				file << setw(15) << setprecision(3) << element_aperture;
				file << setw(10) << shear_failure   << endl;
				*/
			}
		}
	}
	file.close();
	TextInterface fractureAperture;
	fractureAperture.OutputDataAsTextColumns("fractures", *model_, (modelName_ + "_aperture").c_str(), "volume modifier");
}
