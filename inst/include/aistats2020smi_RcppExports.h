// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_aistats2020smi_RCPPEXPORTS_H_GEN_
#define RCPP_aistats2020smi_RCPPEXPORTS_H_GEN_

#include <RcppArmadillo.h>
#include <Rcpp.h>

namespace aistats2020smi {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("aistats2020smi", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("aistats2020smi", "_aistats2020smi_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in aistats2020smi");
            }
        }
    }

    inline Rcpp::List SMI_post_biased_data(const arma::vec Z, const arma::vec Y, const double sigma_z, const double sigma_y, const double sigma_phi, const double sigma_theta, const double sigma_theta_tilde, const double eta = 1) {
        typedef SEXP(*Ptr_SMI_post_biased_data)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_SMI_post_biased_data p_SMI_post_biased_data = NULL;
        if (p_SMI_post_biased_data == NULL) {
            validateSignature("Rcpp::List(*SMI_post_biased_data)(const arma::vec,const arma::vec,const double,const double,const double,const double,const double,const double)");
            p_SMI_post_biased_data = (Ptr_SMI_post_biased_data)R_GetCCallable("aistats2020smi", "_aistats2020smi_SMI_post_biased_data");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_SMI_post_biased_data(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(sigma_z)), Shield<SEXP>(Rcpp::wrap(sigma_y)), Shield<SEXP>(Rcpp::wrap(sigma_phi)), Shield<SEXP>(Rcpp::wrap(sigma_theta)), Shield<SEXP>(Rcpp::wrap(sigma_theta_tilde)), Shield<SEXP>(Rcpp::wrap(eta)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List SMI_pred_biased_data(const arma::vec Z, const arma::vec Y, const double sigma_z, const double sigma_y, const double sigma_phi, const double sigma_theta, const double sigma_theta_tilde, const double eta = 1) {
        typedef SEXP(*Ptr_SMI_pred_biased_data)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_SMI_pred_biased_data p_SMI_pred_biased_data = NULL;
        if (p_SMI_pred_biased_data == NULL) {
            validateSignature("Rcpp::List(*SMI_pred_biased_data)(const arma::vec,const arma::vec,const double,const double,const double,const double,const double,const double)");
            p_SMI_pred_biased_data = (Ptr_SMI_pred_biased_data)R_GetCCallable("aistats2020smi", "_aistats2020smi_SMI_pred_biased_data");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_SMI_pred_biased_data(Shield<SEXP>(Rcpp::wrap(Z)), Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(sigma_z)), Shield<SEXP>(Rcpp::wrap(sigma_y)), Shield<SEXP>(Rcpp::wrap(sigma_phi)), Shield<SEXP>(Rcpp::wrap(sigma_theta)), Shield<SEXP>(Rcpp::wrap(sigma_theta_tilde)), Shield<SEXP>(Rcpp::wrap(eta)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool logd = false) {
        typedef SEXP(*Ptr_dmvnorm_arma)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_dmvnorm_arma p_dmvnorm_arma = NULL;
        if (p_dmvnorm_arma == NULL) {
            validateSignature("arma::vec(*dmvnorm_arma)(arma::mat,arma::rowvec,arma::mat,bool)");
            p_dmvnorm_arma = (Ptr_dmvnorm_arma)R_GetCCallable("aistats2020smi", "_aistats2020smi_dmvnorm_arma");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_dmvnorm_arma(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(mean)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(logd)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::vec >(rcpp_result_gen);
    }

    inline arma::mat categ_to_dummie(const arma::vec x, const arma::vec categ_x) {
        typedef SEXP(*Ptr_categ_to_dummie)(SEXP,SEXP);
        static Ptr_categ_to_dummie p_categ_to_dummie = NULL;
        if (p_categ_to_dummie == NULL) {
            validateSignature("arma::mat(*categ_to_dummie)(const arma::vec,const arma::vec)");
            p_categ_to_dummie = (Ptr_categ_to_dummie)R_GetCCallable("aistats2020smi", "_aistats2020smi_categ_to_dummie");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_categ_to_dummie(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(categ_x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::mat >(rcpp_result_gen);
    }

    inline double sign(const double val) {
        typedef SEXP(*Ptr_sign)(SEXP);
        static Ptr_sign p_sign = NULL;
        if (p_sign == NULL) {
            validateSignature("double(*sign)(const double)");
            p_sign = (Ptr_sign)R_GetCCallable("aistats2020smi", "_aistats2020smi_sign");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_sign(Shield<SEXP>(Rcpp::wrap(val)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double bounce_limit(double x, const double a, const double b) {
        typedef SEXP(*Ptr_bounce_limit)(SEXP,SEXP,SEXP);
        static Ptr_bounce_limit p_bounce_limit = NULL;
        if (p_bounce_limit == NULL) {
            validateSignature("double(*bounce_limit)(double,const double,const double)");
            p_bounce_limit = (Ptr_bounce_limit)R_GetCCallable("aistats2020smi", "_aistats2020smi_bounce_limit");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_bounce_limit(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(a)), Shield<SEXP>(Rcpp::wrap(b)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double dinvgamma(const double x, const double alpha, const double beta, const unsigned int lg = 0) {
        typedef SEXP(*Ptr_dinvgamma)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_dinvgamma p_dinvgamma = NULL;
        if (p_dinvgamma == NULL) {
            validateSignature("double(*dinvgamma)(const double,const double,const double,const unsigned int)");
            p_dinvgamma = (Ptr_dinvgamma)R_GetCCallable("aistats2020smi", "_aistats2020smi_dinvgamma");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_dinvgamma(Shield<SEXP>(Rcpp::wrap(x)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(lg)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline Rcpp::List mcmc_PO(const arma::colvec Y, const arma::mat X, arma::mat X_eta, const unsigned int K, const double power_w, const std::string prior_spec, const bool rnd_eff, const unsigned int n_iter, arma::colvec theta, const arma::mat theta_min_max, arma::colvec theta_prop_int, const std::string theta_prop_kernel, const bool keep_ll = true, const bool check_mcmc = false, const bool verbose = false, const bool quiet = false) {
        typedef SEXP(*Ptr_mcmc_PO)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mcmc_PO p_mcmc_PO = NULL;
        if (p_mcmc_PO == NULL) {
            validateSignature("Rcpp::List(*mcmc_PO)(const arma::colvec,const arma::mat,arma::mat,const unsigned int,const double,const std::string,const bool,const unsigned int,arma::colvec,const arma::mat,arma::colvec,const std::string,const bool,const bool,const bool,const bool)");
            p_mcmc_PO = (Ptr_mcmc_PO)R_GetCCallable("aistats2020smi", "_aistats2020smi_mcmc_PO");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mcmc_PO(Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(X_eta)), Shield<SEXP>(Rcpp::wrap(K)), Shield<SEXP>(Rcpp::wrap(power_w)), Shield<SEXP>(Rcpp::wrap(prior_spec)), Shield<SEXP>(Rcpp::wrap(rnd_eff)), Shield<SEXP>(Rcpp::wrap(n_iter)), Shield<SEXP>(Rcpp::wrap(theta)), Shield<SEXP>(Rcpp::wrap(theta_min_max)), Shield<SEXP>(Rcpp::wrap(theta_prop_int)), Shield<SEXP>(Rcpp::wrap(theta_prop_kernel)), Shield<SEXP>(Rcpp::wrap(keep_ll)), Shield<SEXP>(Rcpp::wrap(check_mcmc)), Shield<SEXP>(Rcpp::wrap(verbose)), Shield<SEXP>(Rcpp::wrap(quiet)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline Rcpp::List mcmc_PO_HM_powered(const arma::mat data_arc, const arma::mat data_mod, double power_w_PO, double power_w_HM, const std::string prior_spec_PO, const std::string prior_spec_HM, const bool PO_site_rnd_eff, const bool HM_site_rnd_eff, const unsigned int n_iter, const unsigned int n_warmup, const unsigned int n_thin, arma::colvec theta, const arma::mat theta_min_max, arma::colvec theta_prop_int, const std::string theta_prop_kernel, arma::colvec ManureLevel_imp, arma::colvec Rainfall_imp, const bool imp_playpen = false, const bool gibbs_hm = true, const bool keep_imp = true, const bool keep_ll = true, const bool check_mcmc = false, const bool verbose = false) {
        typedef SEXP(*Ptr_mcmc_PO_HM_powered)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_mcmc_PO_HM_powered p_mcmc_PO_HM_powered = NULL;
        if (p_mcmc_PO_HM_powered == NULL) {
            validateSignature("Rcpp::List(*mcmc_PO_HM_powered)(const arma::mat,const arma::mat,double,double,const std::string,const std::string,const bool,const bool,const unsigned int,const unsigned int,const unsigned int,arma::colvec,const arma::mat,arma::colvec,const std::string,arma::colvec,arma::colvec,const bool,const bool,const bool,const bool,const bool,const bool)");
            p_mcmc_PO_HM_powered = (Ptr_mcmc_PO_HM_powered)R_GetCCallable("aistats2020smi", "_aistats2020smi_mcmc_PO_HM_powered");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_mcmc_PO_HM_powered(Shield<SEXP>(Rcpp::wrap(data_arc)), Shield<SEXP>(Rcpp::wrap(data_mod)), Shield<SEXP>(Rcpp::wrap(power_w_PO)), Shield<SEXP>(Rcpp::wrap(power_w_HM)), Shield<SEXP>(Rcpp::wrap(prior_spec_PO)), Shield<SEXP>(Rcpp::wrap(prior_spec_HM)), Shield<SEXP>(Rcpp::wrap(PO_site_rnd_eff)), Shield<SEXP>(Rcpp::wrap(HM_site_rnd_eff)), Shield<SEXP>(Rcpp::wrap(n_iter)), Shield<SEXP>(Rcpp::wrap(n_warmup)), Shield<SEXP>(Rcpp::wrap(n_thin)), Shield<SEXP>(Rcpp::wrap(theta)), Shield<SEXP>(Rcpp::wrap(theta_min_max)), Shield<SEXP>(Rcpp::wrap(theta_prop_int)), Shield<SEXP>(Rcpp::wrap(theta_prop_kernel)), Shield<SEXP>(Rcpp::wrap(ManureLevel_imp)), Shield<SEXP>(Rcpp::wrap(Rainfall_imp)), Shield<SEXP>(Rcpp::wrap(imp_playpen)), Shield<SEXP>(Rcpp::wrap(gibbs_hm)), Shield<SEXP>(Rcpp::wrap(keep_imp)), Shield<SEXP>(Rcpp::wrap(keep_ll)), Shield<SEXP>(Rcpp::wrap(check_mcmc)), Shield<SEXP>(Rcpp::wrap(verbose)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<Rcpp::List >(rcpp_result_gen);
    }

    inline arma::colvec loglik_PO_i_cpp(const arma::colvec Y, const arma::mat X, const arma::colvec alpha, const arma::colvec beta) {
        typedef SEXP(*Ptr_loglik_PO_i_cpp)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_loglik_PO_i_cpp p_loglik_PO_i_cpp = NULL;
        if (p_loglik_PO_i_cpp == NULL) {
            validateSignature("arma::colvec(*loglik_PO_i_cpp)(const arma::colvec,const arma::mat,const arma::colvec,const arma::colvec)");
            p_loglik_PO_i_cpp = (Ptr_loglik_PO_i_cpp)R_GetCCallable("aistats2020smi", "_aistats2020smi_loglik_PO_i_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_loglik_PO_i_cpp(Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::colvec >(rcpp_result_gen);
    }

    inline double loglik_PO_cpp(const arma::colvec Y, const arma::mat X, const arma::colvec alpha, const arma::colvec beta) {
        typedef SEXP(*Ptr_loglik_PO_cpp)(SEXP,SEXP,SEXP,SEXP);
        static Ptr_loglik_PO_cpp p_loglik_PO_cpp = NULL;
        if (p_loglik_PO_cpp == NULL) {
            validateSignature("double(*loglik_PO_cpp)(const arma::colvec,const arma::mat,const arma::colvec,const arma::colvec)");
            p_loglik_PO_cpp = (Ptr_loglik_PO_cpp)R_GetCCallable("aistats2020smi", "_aistats2020smi_loglik_PO_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_loglik_PO_cpp(Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double logprior_PO_cpp(const arma::colvec alpha, const arma::colvec beta, const double sigma_eta, const arma::colvec eta, const std::string prior_spec) {
        typedef SEXP(*Ptr_logprior_PO_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_logprior_PO_cpp p_logprior_PO_cpp = NULL;
        if (p_logprior_PO_cpp == NULL) {
            validateSignature("double(*logprior_PO_cpp)(const arma::colvec,const arma::colvec,const double,const arma::colvec,const std::string)");
            p_logprior_PO_cpp = (Ptr_logprior_PO_cpp)R_GetCCallable("aistats2020smi", "_aistats2020smi_logprior_PO_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_logprior_PO_cpp(Shield<SEXP>(Rcpp::wrap(alpha)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(sigma_eta)), Shield<SEXP>(Rcpp::wrap(eta)), Shield<SEXP>(Rcpp::wrap(prior_spec)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline arma::colvec loglik_HM_i_cpp(const arma::colvec Y, const arma::mat X, const arma::colvec beta, const double sigma, const double v, const arma::colvec ind_v) {
        typedef SEXP(*Ptr_loglik_HM_i_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_loglik_HM_i_cpp p_loglik_HM_i_cpp = NULL;
        if (p_loglik_HM_i_cpp == NULL) {
            validateSignature("arma::colvec(*loglik_HM_i_cpp)(const arma::colvec,const arma::mat,const arma::colvec,const double,const double,const arma::colvec)");
            p_loglik_HM_i_cpp = (Ptr_loglik_HM_i_cpp)R_GetCCallable("aistats2020smi", "_aistats2020smi_loglik_HM_i_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_loglik_HM_i_cpp(Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(v)), Shield<SEXP>(Rcpp::wrap(ind_v)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<arma::colvec >(rcpp_result_gen);
    }

    inline double loglik_HM_cpp(const arma::colvec Y, const arma::mat X, const arma::colvec beta, const double sigma, const double v, const arma::colvec ind_v) {
        typedef SEXP(*Ptr_loglik_HM_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_loglik_HM_cpp p_loglik_HM_cpp = NULL;
        if (p_loglik_HM_cpp == NULL) {
            validateSignature("double(*loglik_HM_cpp)(const arma::colvec,const arma::mat,const arma::colvec,const double,const double,const arma::colvec)");
            p_loglik_HM_cpp = (Ptr_loglik_HM_cpp)R_GetCCallable("aistats2020smi", "_aistats2020smi_loglik_HM_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_loglik_HM_cpp(Shield<SEXP>(Rcpp::wrap(Y)), Shield<SEXP>(Rcpp::wrap(X)), Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(v)), Shield<SEXP>(Rcpp::wrap(ind_v)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

    inline double logprior_HM_cpp(const arma::colvec beta, const double sigma, const double v, const double sigma_eta, const arma::colvec eta, const std::string prior_spec) {
        typedef SEXP(*Ptr_logprior_HM_cpp)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP);
        static Ptr_logprior_HM_cpp p_logprior_HM_cpp = NULL;
        if (p_logprior_HM_cpp == NULL) {
            validateSignature("double(*logprior_HM_cpp)(const arma::colvec,const double,const double,const double,const arma::colvec,const std::string)");
            p_logprior_HM_cpp = (Ptr_logprior_HM_cpp)R_GetCCallable("aistats2020smi", "_aistats2020smi_logprior_HM_cpp");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_logprior_HM_cpp(Shield<SEXP>(Rcpp::wrap(beta)), Shield<SEXP>(Rcpp::wrap(sigma)), Shield<SEXP>(Rcpp::wrap(v)), Shield<SEXP>(Rcpp::wrap(sigma_eta)), Shield<SEXP>(Rcpp::wrap(eta)), Shield<SEXP>(Rcpp::wrap(prior_spec)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (Rcpp::internal::isLongjumpSentinel(rcpp_result_gen))
            throw Rcpp::LongjumpException(rcpp_result_gen);
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(Rcpp::as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<double >(rcpp_result_gen);
    }

}

#endif // RCPP_aistats2020smi_RCPPEXPORTS_H_GEN_