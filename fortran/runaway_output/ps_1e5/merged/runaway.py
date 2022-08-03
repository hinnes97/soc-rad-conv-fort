import numpy as np
import phys
from moistad import calc_q, moistadiabat, satvp, dew_point, cold_trap, dlntdlnp
import scipy.integrate as spi
import scipy.optimize as spo
import scipy.interpolate as interp
import matplotlib.pyplot as plt
from soc_OLR import calc_OLR_soc, calc_SW_soc
import params
import sys

def calc_tau(p, T, q):
    """Calculate tau as given in Ding 2016. Uses constant kappa for water
       dependent on water path, plus a constant background term tau_0 for 
       the dry gasses"""

    # Constant term
    tau_dry = (p/p[-1])*params.tau_0

    #q = calc_q(T, p, params.dry, params.wet)

    
    tau_wet = np.zeros_like(p)
    tau_dry = np.zeros_like(p)
    tau_dry_1 = np.zeros_like(p)
    q_mid = np.sqrt(q[1:]*q[:-1])
    p_mid = (p[1:]+ p[:-1])/2
    dp = np.diff(p)

    R = q*params.Rc + (1-q)*params.dry.R
    R_mid= (R[1:] + R[:-1])/2
    T_mid = (T[1:] + T[:-1])/2
    for i in range(1,len(p)):
        tau_wet[i] = tau_wet[i-1] + dp[i-1]*params.kappa/params.grav*q_mid[i-1]*2
        tau_dry_1[i] = tau_dry_1[i-1] + dp[i-1]*(1-q_mid[i-1])*2/params.grav*params.kap0
        #print(dp[i-1],(1-q_mid[i-1]),params.grav,params.kap0)
        tau_dry[i] = tau_dry[i-1] + p_mid[i-1]/R_mid[i-1]/T_mid[i-1] * params.k2 * dp[i-1]/params.grav*(1-q_mid[i-1])**2 

#    print('linear:', tau_dry_1[-1], 'CIA', tau_dry[-1])
    #print((tau_wet/tau_dry)[-1])
    
    tau = tau_dry + tau_wet + tau_dry_1

    #i = min(np.searchsorted(tau, 1), len(tau)-1)
    #print(tau_dry_1[i], tau_wet[i], tau_wet[i]/tau_dry_1[i], tau_dry[i])
    #print(tau[-1])
    return tau_wet + tau_dry_1 + tau_dry


def calc_OLR_toon(p, T, q):
    """Uses Toon method of calculating flux, but simplified for pure radiative atmosphere
       Equation taken from Heng 2014 (eq. 47). Switches to the diffusion approximation if 
       tau >> 1"""

    F0 = phys.sigma*T[-1]**4

    tau = calc_tau(p,T,q)
    
    # Don't bother integrating when tau>100, here use diffusion approximation
    # In reality only need the boundary condition, since flux just depends on local dT/dp
    taulim = 100
    dTdp = T/p * np.gradient(np.log(T), np.log(p))

    if tau[-1]<taulim:
        T_rad = T
        P_rad = p
        tau_rad = tau
        flup_boundary = phys.sigma*T[-1]**4
    else:
        flup_boundary = phys.sigma * T[tau<taulim][-1]**4 

        # Use this for toon method (use characteristic slant path mu = 1/2)
        T_rad = T[tau<taulim]
        P_rad = p[tau<taulim]
        tau_rad = tau[tau<taulim]
    mu = 1/2
    B_dash = phys.sigma*(T_rad[:-1]**4 - T_rad[1:]**4)/mu/np.diff(tau_rad)
    B1_plus = phys.sigma*T_rad[1:]**4 + B_dash*mu*(-np.diff(tau_rad)) + B_dash/2.
    B2_plus = phys.sigma*T_rad[1:]**4 + B_dash/2.

    trans = np.exp(-(2*mu)*np.diff(tau_rad))

    olr = flup_boundary
    for k in range(len(T_rad)-2, -1, -1):
        olr = olr*trans[k] + B1_plus[k] - B2_plus[k]*trans[k]

    contrib = None

    #print('sig*T^4 at tau=1', phys.sigma*T[np.argmin(np.absolute(tau-1))]**4)
    #print('OLR: ', olr)
    return olr, contrib


    
def crit_limit(t, y, S, pvec, qmin):
    # Check if under the critical temperature
    q = calc_q(np.exp(y[0]), np.exp(t), params.dry, params.wet)
    q = min(qmin, q)

    q_crit = np.exp(y[0])*params.Rc/params.L/params.pi
    #print(t,y,q, q_crit)
    
    return (q - q_crit)

def adrad_limit(t, y, S, pvec, qmin):
    ad_grad = dlntdlnp(t, y, params.dry, params.wet)
    rad_grad = new_rad_gradient(t,y,S,pvec,qmin) 

    return rad_grad - ad_grad

def calc_OLR(p, T, q):
    F0 = phys.sigma*T[-1]**4 # Surface term

    # Calculate tau
    tau = calc_tau(p, T, q)
    
    # Calculate integrated part
    #y = phys.sigma*T**4*np.exp(-tau)
    #int_term = spi.trapezoid(y, tau)
    y = -phys.sigma*T**4
    int_term = spi.trapezoid(y, np.exp(-tau))
        
    #np.savetxt('rad.txt', [[tt,yy] for (tt,yy) in zip(tau,y)])
    OLR = F0*np.exp(-tau[-1]) + int_term
    #('int term', int_term, 'bound term', F0*np.exp(-tau[-1]))
    D=1
    contrib = phys.sigma*(T[1:]/2+T[:-1]/2)**4*np.pi*2/D *np.diff(np.exp(-D*(tau[-1] - tau)))
    #print('DTAU', np.diff(tau))
    return OLR, contrib

def diffusion_gradient(p, T, S, pvec, qmin):
    """Calculates the gradient in the diffusion approximation limit
       used to calculate the lapse rate when convection is inhibited
       Flux = 16/3 * g/k * dT/dp * sigma*T**3. For atmosphere heated
       at the surface, flux = instellation, S"""


    q = calc_q(T,p, params.dry, params.wet)
    if params.old_scheme:
        pass
    else:
        qmin = min(q, qmin)
        q = qmin
    #if params.old_scheme:
    #    q = calc_q(T, p, params.dry, params.wet)
    #else:
    #    pass
    #kappa_0 = params.tau_0/params.mass_path
    
    index = np.argmin(np.absolute(p-pvec))
    S_1 = S[index]
    R = q*params.Rc + (1-q)*params.dry.R

    dTdp = 3./16. * (params.kappa*q+params.kap0*(1-q) + params.k2*(1-q)**2*(p/R/T))/params.grav * S_1 /phys.sigma/T**3

    #dTdp = 3./16. * (params.kappa*q+params.kap0*(1-q))/params.grav * S_1 /phys.sigma/T**3
    return dTdp


def new_ad_gradient(logp, logT, S, pvec, qmin):
    #p = np.exp(logp)
    #T = np.exp(logT)[0]
    
    ad_grad = dlntdlnp(logp, logT[0], params.dry, params.wet)
    
    return ad_grad

def new_rad_gradient(logp, logT, S, pvec, qmin):

    p = np.exp(logp)
    
    T = np.exp(logT)[0]
    #T  = T[0]
    q = calc_q(T,p, params.dry, params.wet)

    qmin = min(q, qmin)
    q = qmin

    #index = np.argmin(np.absolute(p-pvec))
    #S_1 = S[index]
    #S_1 = S
    S_1 = S*np.exp(-params.kappa_sw*p/params.grav)

    #S_1 = S*np.exp(-params.kappa_sw*p/params.grav)
    R = q*params.Rc + (1-q)*params.dry.R

    rad_grad = 4./16.  * p * (params.kappa*q+params.kap0*(1-q) + params.k2*(1-q)**2*(p/R/T))/params.grav * S_1 /phys.sigma/T**4

    #print('Solar', S, 'rad_grad', rad_grad)
    return rad_grad

def solve_adrad_region(pbot, tbot, qbot, S, p_eval, func, events=(crit_limit, )):

    # Change into log coordinates

    logpbot = np.log(pbot)
    logtbot = np.log(tbot)
    logp_eval = np.log(p_eval)
    #logpbot = pbot
    #logtbot = tbot
    #logp_eval = p_eval
    
    q = qbot
    logt = logtbot
    qmin = q
    
    events, is_terminal, event_dir = spi._ivp.ivp.prepare_events(events)

    if events is not None:
        events = [lambda p,T, event=event : event(p,T,S,p_eval,qmin) for event in events]
        g = [event(logpbot,logtbot) for event in events]
        logp_events = [[] for _ in range(len(events))]
        logt_events = [[] for _ in range(len(events))]
    else:
        logt_events = None
        logp_events = None


    fun = lambda logp,logT,func = func: func(logp,logT,S, p_eval, qmin)
    solver = spi.DOP853(fun, logpbot, logtbot, np.log(params.pup), max_step = 0.1)

    logps = []
    logts = []

    logp_eval_i = p_eval.shape[0]



    #params.rad_region = False
    while solver.status == 'running':
        
        #print('before step')
        message = solver.step()
        #print('after step')
        
        #print('current pressure', solver.t)
        if message == 'finished':
            return np.flip(np.exp(np.array(ts[0]))), np.exp(np.flip(np.array(ps)))
        elif message == 'failed':
            exit('failure in solver')

        sol = None
        
        logp_old = solver.t_old
        logp = solver.t
        logt_old = logt
        logt = solver.y
        terminate = False
        #print(logp, logt[0])
        if events is not None:
            # Check for events
            g_new = [event(logp, logt) for event in events]

            active_events = spi._ivp.ivp.find_active_events(g, g_new, event_dir)

            
            if active_events.size>0:
                
                #print('active_event')
                if sol is None:
                    sol = solver.dense_output()

                # Find the root not using the IVP method


                root_indices, roots, terminate = spi._ivp.ivp.handle_events(
                        sol, events, active_events, is_terminal, logp_old, logp)

                for e, pe in zip(root_indices, roots):
                    logp_events[e].append(pe)
                    logt_events[e].append(sol(pe))

                #print('event activated', np.exp(logp), np.exp(logt))
                if terminate:
                    status = 1
                    logp = roots[-1]
                    logt = sol(logp)

                    # Work out if we should print to output
                    logp_eval_i_new = np.searchsorted(np.log(p_eval), logp, side='left')
                    logp_eval_step = logp_eval[logp_eval_i_new:logp_eval_i][::-1]

                    if logp_eval_step.size >0:
                        sol = solver.dense_output()
                        logps.append(logp_eval_step)
                        #print('p_eval_step', p_eval_step)
                        logts.append(sol(logp_eval_step))


                    #ps.append([p])
                    #ts.append([t])


                    logts = np.hstack(logts)
                    logps = np.hstack(logps)

                    logp_events = [np.asarray(pe) for pe in logp_events]
                    logt_events = [np.asarray(te) for te in logt_events]
                    
#                    print('terminating', terminate, np.exp(logp), np.exp(logps[-1]))
                    return np.flip(np.exp(logts[0])), np.flip(np.exp(logps)), logt_events, logp_events,calc_q(np.exp(logt[0]), np.exp(logp), phys.H2, phys.H2O), terminate
                    #return np.flip(ts[0]), np.flip(ps), logt_events, logp_events, terminate


            g = g_new

        #print('shouldnt end up here')
        # Work out if we should print to output
        logp_eval_i_new = np.searchsorted(np.log(p_eval), logp, side='left')
        logp_eval_step = logp_eval[logp_eval_i_new:logp_eval_i][::-1]

        if logp_eval_step.size >0:
            sol = solver.dense_output()
            logps.append(logp_eval_step)
            #print('p_eval_step', p_eval_step)
            logts.append(sol(logp_eval_step))
            logp_eval_i = logp_eval_i_new
        
        # Work out new q
        q_new = calc_q(np.exp(logt[0]),np.exp(logp),phys.H2, phys.H2O)
        #q_new = calc_q(logt, logp, phys.H2, phys.H2O)
        # Require q to be less than old q (cold trap)
        qmin = min(q_new, qmin)
        q = qmin

    return np.flip(np.exp(np.hstack(logts)[0])), np.flip(np.exp(np.hstack(logps))), logt_events, logp_events, q, terminate




def grey_gas_solution(p, p0, T0):
    # Calculate tau from TOA down to the point p0, T0
    twosigt4 = 2*phys.sigma*T0**4

    tau = calc_tau()
    # Downwelling SW radiation
    #Sm = 
    int_term  = trapezoid()


def mixed_limits(logp, logT, S, pvec, qmin):
    
    return adrad_limit(logp, logT, S, pvec, qmin)*crit_limit(logp, logT, S, pvec, qmin)
    
def find_mixed_profile_streamlined(S, Ts, ps):
    #print('IN MIXED PROFILE STREAMLINED')
    P = np.logspace(np.log10(params.pup), np.log10(ps), params.N)

    q_bot = calc_q(Ts, ps, params.dry, params.wet)

    p_start = P[-1]
    T_start = Ts
    Pevals = P

    psols=[]
    tsols=[]
    crit_limit.terminal = True
    adrad_limit.terminal = True


    if params.pure_radiative == True:
        T,p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_rad_gradient,
                                                                         events = None)
                        
        tsols = np.hstack((T, tsols))
        psols = np.hstack((p, psols))
        
        return psols, tsols, []

        
    events = []
    evnts = (crit_limit, adrad_limit)
    event_type = ''

    qlim = crit_limit(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot)
    adrad = adrad_limit(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot)

    if params.CI == False:
        if adrad > 0:
            stable = False
            evnts = (adrad_limit,)
            rad_to_ad = None
        else:
            stable = True
            evnts = None
    else:
        #print(mixed_limits(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot))
        #print(new_ad_gradient(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot), new_rad_gradient(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot))
        #print('HELLO WORLD', S, S*np.exp(-params.kappa_sw*p_start/params.grav))
        if mixed_limits(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot) > 0:
            stable = True
        else:
            stable = False

    while True:
        if stable:
            fun = new_rad_gradient
        else:
            fun = new_ad_gradient

        if (not params.CI) and (not params.rad_strat):
            # Just pure adiabat here
            evnts = None
            fun = new_ad_gradient
#        elif (not params.CI) and event_type=='gradient':
            # Should only stop when gradient changes
#            evnts = (adrad_limit,)
#            fun = new_ad_gradient
#            adrad_limit.direction = -1
#            rad_to_ad = None
        elif (not params.rad_strat):
            # Make sure if under qc then on adiabat
            if (event_type=='' and qlim < 0) or event_type == 'moisture':
                fun = new_ad_gradient
                evnts = None
            

        T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         fun,
                                                                         events = evnts)
                
        tsols = np.hstack((T, tsols))
        psols = np.hstack((p, psols))

        if terminal:
            if len(p_events)>1:
                if len(p_events[0])>0:
                    i = 0
                    event_type = 'moisture'
                    evnts = (adrad_limit,)
                    rad_to_ad = None
                else:
                    i = 1
                    event_type= 'gradient'
                    evnts = (crit_limit,)
                    rad_to_ad = stable

            else:
                if rad_to_ad is None:
                    event_type = 'gradient'

                    if params.CI:
                        evnts = (crit_limit,)
                    else:
                        evnts = None
                    rad_to_ad = stable
                    i=0

                else:
                    event_type = 'moisture'
                    evnts = (adrad_limit,)
                    rad_to_ad = None
                    i=0



            event_dict = {'type': event_type,
                          'pres': np.exp(p_events[i][0]),
                          'temp': np.exp(t_events[i][0][0]),
                          'q'   : q_event,
                          'rad_to_ad':rad_to_ad}
            events.append(event_dict)

            stable = not stable
            T_start = np.exp(t_events[i][0][0])
            p_start = np.exp(p_events[i][0])
            q_bot = q_event

            Pevals = np.r_[P[P<p_start], p_start]
        else:
            return psols, tsols, events
                

def find_mixed_profile_new(S, Ts, ps):
    P = np.logspace(np.log10(params.pup), np.log10(ps), params.N)

    q_bot = calc_q(Ts, ps, params.dry, params.wet)
    crit_limit.terminal = True
    adrad_limit.terminal = True
    
    # Check crit limit function to see what to start with
    q_check = crit_limit(np.log(ps), [np.log(Ts)], S, P, q_bot)
    grad_check = adrad_limit(np.log(ps), [np.log(Ts)], S, P, q_bot)

    p_start = P[-1]
    T_start = Ts
    Pevals = P

    psols = []
    tsols = []

    #print(q_check, grad_check)
    #print('beginning mixed profile', p_start)
    if q_check>0:
        if grad_check > 0:
            #print('q>0, gradcheck>0')
            #T, p, t_events, p_events,  terminal = solve_rad_region(p_start, [T_start],
            #                                                             q_bot, S, Pevals,
            #                                                             events = (crit_limit, adrad_limit))
            T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_rad_gradient,
                                                                         events = (crit_limit, adrad_limit))

            tsols = np.hstack((T, tsols))
            psols = np.hstack((p, psols))

            if terminal==True:
                # Find out what kind of event it was
                if len(p_events[0])>0:
                    # reached critical q
                    T_start = np.exp(t_events[0][0][0])
                    p_start = np.exp(p_events[0][0])
                    q_bot = q_event

                    Pevals =np.r_[P[P<p_start], p_start]
                    grad_check = 1
                    #print('Switching from q>qc rad to q<qc ad')
                    #print('p', p_start, 't', T_start,'q', q_bot)

                else:
                    T_start = np.exp(t_events[1][0][0])
                    p_start = np.exp(p_events[1][0])
                    q_bot = q_event

                    Pevals = np.r_[P[P<p_start], p_start]

                    #print('Switching from q>qc rad to q>qc ad')
                    #print('p', p_start, 't', T_start,'q', q_bot)
                    T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_ad_gradient,
                                                                         events = (crit_limit))
                    tsols = np.hstack((T, tsols))
                    psols = np.hstack((p, psols))

                    T_start = np.exp(t_events[0][0][0])
                    p_start = np.exp(p_events[0][0])
                    Pevals = np.r_[P[P<p_start], p_start]
                    q_bot = q_event
                    grad_check = -1
                    #print('Switching from q>qc ad to q<qc rad')
                    #print('p_start', p_start, 'T_start', T_start, 'psols[0]', psols[0])
            else:
                return psols, tsols
        else:
            #print('q>0, grad_check<0')
            
            T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_ad_gradient,
                                                                         events = (crit_limit))

            tsols = np.hstack((T, tsols))
            psols = np.hstack((p, psols))
            
            if terminal == True:
                T_start = np.exp(t_events[0][0][0])
                p_start = np.exp(p_events[0][0])
                
                q_bot = q_event

                Pevals = np.r_[P[P<p_start], p_start]

                grad_check = -1
                #print('Switching from q>qc ad to q<qc rad')
                #print('p', p_start, 't', T_start,'q', q_bot)
                #print('qstart', calc_q(T_start, p_start, params.dry, params.wet), params.Rc*T_start/params.L/params.pi)
            else:
                return psols, tsols


    # Now we are below crit_limit q<qc
    #print('grad_check check', grad_check)
    if grad_check < 0:
        #T, p, t_events, p_events, terminal = solve_rad_region(p_start, [T_start],
        #                                                                 q_bot, S, Pevals,
        #                                                                 events = None)
        T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_rad_gradient,
                                                                         events = (None))

        tsols = np.hstack((T, tsols))
        psols = np.hstack((p, psols))
    else:
        #print('grad_check should be >0 ')
        #print('p_start', p_start, 'T_start', T_start)
        #print('grad_check', adrad_limit(np.log(p_start), [np.log(T_start)], S, Pevals, q_bot))
        T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_ad_gradient,
                                                                         events = adrad_limit)
        tsols = np.hstack((T, tsols))
        psols = np.hstack((p, psols))

        
        if terminal==True:

            T_start = np.exp(t_events[0][0][0])
            p_start = np.exp(p_events[0][0])
            Pevals = np.r_[P[P<p_start], p_start]
            q_bot = q_event
            #print('Switching from q<qc ad to q<qc rad')
            #print('p', p_start, 't', T_start,'q', q_bot, 'Pevals', Pevals)

            T, p, t_events, p_events, q_event, terminal = solve_adrad_region(p_start, [T_start],
                                                                         q_bot, S, Pevals,
                                                                         new_rad_gradient,
                                                                             events = None)

            #T, p, t_events, p_events,  terminal = solve_rad_region(p_start, [T_start],
            #                                                             q_bot, S, Pevals,
            #                                                             events = None)

            tsols = np.hstack((T, tsols))
            psols = np.hstack((p, psols))

    #print('ending new profile')
    return psols, tsols

    # while True:
    #     if check>0:
    #         func = new_rad_gradient
    #         crit_limit.direction = -1
    #     else:
    #         func = new_ad_gradient
    #         crit_limit.direction = 1

    #     T, p, t_events, p_events, q_event,  terminal = solve_adrad_region(p_start, [T_start], qbot, S, Pevals,func,  events=(crit_limit,))

    #     psols = np.hstack((p, psols))
    #     tsols = np.hstack((T, tsols))
    #     qbot = q_event - q_event*0.00001
        
    #     if terminal==False:
    #         #print(qbot, params.Rc*Ts/params.L/params.pi, params.Rc*tsols[0]/params.L/params.pi)
    #         break
    #     else:
    #         T_start = t_events[0][0][0] 
    #         p_start = p_events[0][0] 
    #         Pevals = np.r_[P[P<p_start], p_start]
    #         print('switching', qbot, T_start*params.Rc/params.L/params.pi, T_start, p_start)
    #         print(new_ad_gradient(np.log(p_start), np.log(T_start),  S, Pevals,q_event), new_rad_gradient(np.log(p_start), np.log(T_start), S, Pevals, q_event))
    #         check = check*-1
    # return psols, tsols
    
    
def OLR_diff(S, Ts, ps):

    # Find S in the radiative layer using simple exponential decay
    #P_tot, T_tot, p_strat = find_mixed_profile(S, Ts, ps)
    #print('OLR DIFF', S)
    P_tot, T_tot, events = find_mixed_profile_streamlined(S, Ts, ps)
   
    T_tot, q = cold_trap(T_tot, P_tot, params.dry, params.wet, params.strat_temp)

    if params.soc:
        OLR, y = calc_OLR_soc(P_tot, T_tot, q)
    else:
        OLR,y = calc_OLR(P_tot, T_tot, q)
#    print('S', S, 'OLR', OLR)
    return S - OLR


def three_layer_atm(p_steam, Ts):
    """Pure steam layer up to p_steam, then radiative layer until q_crit, then moist
       adiabatic layer on top of that"""
    # Pure steam layer
    ps = satvp(Ts, params.wet)

    #print('PS, ', ps, 'P_steam', p_steam)
    M = 100
    p_steams = np.logspace(np.log10(p_steam), np.log10(ps), M)
    T_steams = dew_point(p_steams, params.wet)
    #print(p_steams, T_steams)

    # Call find_mixed_profile, using p_steams, T_steams as BC
    #print('p_steam', p_steam, 'T_steams[0]', T_steams[0])
    print('low S', OLR_diff(100, T_steams[0], p_steam), 'high_s', OLR_diff(1000, T_steams[0], p_steam))

    result = spo.root_scalar(OLR_diff, args=(T_steams[0], p_steam), bracket=[1,1000], rtol=0.00001)
    OLR = result.root
    #print('three layer atm', OLR)
        #S_rad = OLR*np.exp(-params.kappa_sw*params.mass_path)
#        print('S_rad', S_rad)
        #new_P, T = find_mixed_profile(S_rad, Ts, ps)

        #P, T, p_strat = find_mixed_profile(OLR, T_steams[0], p_steam)
    P, T,events = find_mixed_profile_streamlined(OLR, T_steams[0], p_steam)

    P = np.r_[P, p_steams[1:]]
    T = np.r_[T, T_steams[1:]]
    #print(T)
    return P, T, events

def get_profile_only(Ts, ps, three_layer=False, p_steam=1.e5):
    """ Same as get_profile_and_OLR, except don't calculate OLR so often since, 
        when using SOCRATES, this is too expensive a calculation"""

    P = np.logspace(np.log10(params.pup),np.log10(ps), params.N)

    q_crit_s = params.Rc*Ts/params.L/params.pi

    qs = calc_q(Ts, ps, params.dry, params.wet)

    #print('Critical surface q:', q_crit_s)
    #print('Actual surface q:', qs)
    #print('Ts, Ps', Ts, ps)
    # Check Leconte criterion

    if three_layer:
        P, T, events = three_layer_atm(p_steam, Ts)

        T, q = cold_trap(T, P, params.dry, params.wet, params.strat_temp)

        return P, T,q, events
    else:
#        print('COMP INHIBITED')
        # Require root finding s.t. S = OLR, need S for radiative part of atm
    
        if (params.CI == False) and (params.rad_strat ==False):
            Tad,Pad = moistadiabat(P, Ts, params.dry, params.wet)
            Tad, q = cold_trap(Tad, P, params.dry, params.wet, params.strat_temp)
            return P, Tad, q, []
        else:
            result = spo.root_scalar(OLR_diff, args=(Ts, ps), bracket=[1e-8,1000], rtol=0.00001)
            OLR = result.root
            #S_rad = OLR*np.exp(-params.kappa_sw*params.mass_path)

            #new_P, T,  p_strat = find_mixed_profile(OLR, Ts, ps)
            new_P, T,events = find_mixed_profile_streamlined(OLR, Ts, ps)

            T, q = cold_trap(T, new_P, params.dry, params.wet, params.strat_temp)

            return new_P, T, q, events

def get_profile_and_OLR(Ts, ps, three_layer=False, p_steam=1.e5):

    P = np.logspace(np.log10(params.pup),np.log10(ps), params.N)

    q_crit_s = params.Rc*Ts/params.L/params.pi

    qs = calc_q(Ts, ps, params.dry, params.wet)

    #print('Critical surface q:', q_crit_s)
    #print('Actual surface q:', qs)
    #print('Ts, Ps', Ts, ps)
    # Check Leconte criterion
    #print('Ps', ps, 'Ts', Ts, 'qs', qs, 'crit qs', q_crit_s)

    if three_layer:
        P, T, events = three_layer_atm(p_steam, Ts)
        
                
        T, q = cold_trap(T, P, params.dry, params.wet, params.strat_temp)
        
        if params.soc:
            OLR, contrib = calc_OLR_soc(P,T,q)
        else:
            OLR, contrib = calc_OLR(P, T, q)
        
        return OLR, P, T, contrib, q, events
    else:
        #print('COMP INHIB', qs, q_crit_s, ps)
        # Require root finding s.t. S = OLR, need S for radiative part of atm
        result = spo.root_scalar(OLR_diff, args=(Ts, ps), bracket=[1e-8,1000], rtol=0.00001)
        OLR = result.root
        #S_rad = OLR*np.exp(-params.kappa_sw*params.mass_path)
        #new_P, T, p_strat = find_mixed_profile(OLR, Ts, ps)
        new_P, T, events = find_mixed_profile_streamlined(OLR, Ts, ps)
        
        T, q = cold_trap(T, new_P, params.dry, params.wet, params.strat_temp)

        if params.soc:
            OLR, contrib = calc_OLR_soc(new_P,T,q)
        else:
            OLR, contrib = calc_OLR(new_P,T, q)

        return OLR, new_P, T, contrib, q, events

def dry_path(P, T, q):
    path = (1-q)/params.grav
    
    return spi.trapezoid(path, P)

def surface_pressure_finder(log10ps, Ts, mass_path):
    """For a given Ts, mass path, find the ps of the curve that will 
       give the correct mass path"""

    #print('SURFACE PRESSURE ', 10**log10ps)
    P, T, q,events = get_profile_only(Ts, 10**log10ps)
    #print('in root finding procedure', P[100], T[100], q[100])
    m_path = dry_path(P,T, q)
    #print(10**log10ps, mass_path, m_path, q)
    #print('PS', 10**log10ps, 'Ts',Ts, 'MASS PATH', m_path)
    
    #print(mass_path, m_path, mass_path-m_path)
    return mass_path - m_path

def p_steam_finder(p_steam, Ts, mass_path):
    P, T, q,events = get_profile_only(Ts, satvp(Ts, params.wet), True, p_steam)

    m_path = dry_path(P,T, q)

    return mass_path - m_path

def pure_steam_curve(Ts_vector):
    OLRs = np.zeros_like(Ts_vector)
    ps_vec = np.zeros_like(Ts_vector)
    params.tau_0 = 0
    for i,Ts in enumerate(Ts_vector):
        p_sat = satvp(Ts, params.wet)
        ps_vec[i] = p_sat

        P = np.logspace(np.log10(params.pup), np.log10(p_sat), params.N)
        T = dew_point(P, params.wet)

        if params.soc:
            OLRs[i], contrib = calc_OLR_soc(P, T, np.ones_like(T))
            sw_down = calc_SW_soc(P, T, np.ones_like(T))

        else:
            test = calc_OLR(P,T, np.ones_like(T,dtype=float))
            OLRs[i], contrib = calc_OLR(P, T, np.ones_like(T))

        print(Ts, OLRs[i])

    return OLRs, ps_vec

def Ts_vs_OLR(Ts_vector):
    OLRs = np.zeros_like(Ts_vector)
    ps_vec = np.zeros_like(Ts_vector)
    p_steams = np.zeros_like(Ts_vector)
    dry_path_vec = np.zeros_like(Ts_vector)
    events_vec = []
    #p_strats = np.zeros_like(Ts_vector)
    three_layer_vec = np.zeros_like(Ts_vector, dtype=bool)

    do_try_statement = True # Once p_steam is found, it will be identical for rest of Ts vectors, so
                            # don't bother iterating any more
    p_steam = 10*params.pup                        
    for i,Ts in enumerate(Ts_vector):
        p_sat = satvp(Ts, params.wet)
        p_d = params.mass_path*params.grav
        
        # Rough upper limit on surface pressure by solving if q=constant at value of surface
        pmax = 0.5*(p_sat + p_d) + ((p_sat+p_d)**2 + 4*(9.-1.)*p_sat*p_d)**0.5

        #ps = params.mass_path*params.grav + p_sat
        
        if do_try_statement:
            try:
    #            print('MASS PATH SIGN: ', surface_pressure_finder(np.log10(p_sat), Ts, params.mass_path, CI), surface_pressure_finder(np.log10(p_sat), Ts, params.mass_path, CI))
                #print('LOWER', surface_pressure_finder(np.log10(max(p_sat*1.000001, params.mass_path*params.grav)), Ts, params.mass_path), np.log10(max(p_sat*1.000001, params.mass_path*params.grav)))
                #print('UPPER', surface_pressure_finder(8, Ts, params.mass_path))
                 
#                if params.pure_radiative:
#                    raise ValueError
                result = spo.root_scalar(surface_pressure_finder, args=(Ts, params.mass_path), bracket=[np.log10(max(p_sat*1.000001, params.mass_path*params.grav)),8 ], rtol=0.0000001)

                #print(result.iterations, result.function_calls, result.flag)
                ps = 10**result.root
            
                #print('CHECKING SOLUTION', surface_pressure_finder(np.log10(ps), Ts, params.mass_path)) 
                #ps = params.mass_path*params.grav + p_sat
                three_layer = False
                p_steam = ps
                three_layer_vec[i] = three_layer
    
            except ValueError as err:
            #     # Surface pressure with given dry mass does not exist. Put a pure water layer below the radiative layer

                print('THREE LAYER ACTIVATED')
                
                
                ps = satvp(Ts, params.wet)

                try:
                    print(p_steam_finder(1000, Ts, params.mass_path), 'UPPER')
                    print(p_steam_finder(ps, Ts, params.mass_path), 'LOWER')
                    p_steam = spo.root_scalar(p_steam_finder, args=(Ts, params.mass_path),
                                               bracket=[p_steam*0.999, ps*0.999]).root
                except ValueError as e:
                    # Unable to find solution, just return NaNs
                    print(e)
                    print('No further solutions')
                    raise(e)
                    three_layer_vec[i:] = np.nan
                    ps_vec[i:] = np.nan
                    OLRs[i:] = np.nan
                    dry_path_vec[i:] = np.nan
                    p_steams[i:] = np.nan
                    #p_strats[i:] = np.nan
                    return OLRs, ps_vec, dry_path_vec, p_steams, three_layer_vec
                    
                three_layer = True
                do_try_statement=False
                three_layer_vec[i] = three_layer

        else:
            # Once we are above critical temperature, we now have three layers and p_steam
            # remains the same, since atmosphere below p_steam is 100% water and
            # confined to a single steam adiabat 
            
            ps = satvp(Ts, params.wet)
            three_layer_vec[i] = True
                
        OLR, P, T,  contrib, q, events= get_profile_and_OLR(Ts, ps, three_layer, p_steam)
        
        #if params.soc:
            #sw_down = calc_SW_soc(P, T, q)
                
            #print('Qs',np.amax(q), np.amax(params.Rc*T/params.L/params.pi))
                
            #p_crit = (P[q>params.Rc*T/params.L/params.pi])[0]
                
            #print('SW down in rad layer', sw_down[(P<p_steam) & (P>p_crit)])
            #Print('SW down c.f. instellation' ,sw_down[(P<p_steam) & (P>p_crit)]/(1368./4.))
                
            
        print('----------------------------------------------------------------------------------------------------------') 
        print('Ts', Ts,'Ps', ps, 'OLR',  OLR, 'DRY PATH', dry_path(P,T, q))
#        print(OLR, P, T, calc_q(T, P, dry, wet))
        print('----------------------------------------------------------------------------------------------------------') 
        #q = calc_q(T,P,params.dry,params.wet)
        if np.sum(q>1.01)>0:
            print('WARNING, q>1', params.mass_path,'Ts',  Ts, 'PS', ps, (P[q>1])[-1], q[q>1])
    
        ps_vec[i] = ps
        OLRs[i] = OLR
        dry_path_vec[i] = dry_path(P, T, q)
        events_vec.append(events)
        #p_strats[i] = p_strat
        #np.savetxt('dat.txt', [[p, t] for (p,t) in zip(P,T)])
        p_steams[i] = p_steam

    return OLRs, ps_vec, dry_path_vec, p_steams, three_layer_vec, events_vec

def __main__():

    composition_inhibited=True
    Ts_vec = np.linspace(200,500, 100)
    OLR_vec, ps_vec  = Ts_vs_OLR(Ts_vec)

    
    composition_inhibited= False

    OLR, ps = Ts_vs_OLR(Ts_vec, composition_inhibited)
    plt.figure()
    plt.plot(Ts_vec, OLR_vec)
    plt.plot(Ts_vec, OLR)

    ax = plt.twinx()
    ax.semilogy(Ts_vec, ps_vec, color='r')
    ax.semilogy(Ts_vec, ps, color='g')
    
    plt.xlabel('Surface temperature (K)')
    plt.ylabel('OLR (W/m^2)')
    plt.tight_layout()
    plt.savefig('Ts_vs_OLR_2.pdf')

    
#    OLR, P, Tad = get_profile_and_OLR(300, 1.e7, adiabat=True)
#    q = calc_q(Tad, P, dry, wet)
#    OLR_2, P_2, Tad_2 = get_profile_and_OLR(300, 1.e7, adiabat=False)
#    q_2 = calc_q(Tad_2, P_2, dry, wet)

    # print(P, Tad)
    # print(P_2, Tad_2)
    # print('OLR = ', OLR, 'OLR_2 = ', OLR_2)
#    print('DRY PATH 1=  ', dry_path(P, Tad), 'DRY PATH 2= ', dry_path(P_2, Tad_2))
    # plt.figure()
    # #plt.semilogy(Tad, P, Tad_2, P_2)
    # plt.scatter(Tad, P)
    # plt.scatter(Tad_2, P_2)
    # plt.gca().set_yscale('log')
    # plt.gca().invert_yaxis()
    # plt.xlabel('Temperature (K)')
    # plt.ylabel('Pressure (Pa)')
    # plt.tight_layout()
    # plt.savefig('profile_test.pdf')
    
if __name__ == '__main__':
    __main__()
