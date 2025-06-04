using Random
using Statistics

function compute_entropy(alpha; k=3, N=100_000, num_iter=100, num_samples=10_000)
    """
    BP algorithm for calculating entropy density in random 3-SAT
    Parameters:
        alpha : clause density
        k     : SAT problem order (default 3-SAT)
        N     : population size
        num_iter : iterations for population dynamics
        num_samples : Monte Carlo samples for entropy calculation
    """
    # Initialize populations with uniform distributions
    rng = MersenneTwister()  # Separate RNG for thread safety
    eta = rand(rng, N)
    eta_s = rand(rng, N)
    
    # Population dynamics iterations
    for _ in 1:num_iter
        # Update ηs population (clause to variable messages)
        for i in 1:N
            selected = rand(rng, 1:N, k-1)
            new_eta_s = 1.0 - prod(1.0 .- eta[selected])
            replace_idx = rand(rng, 1:N)
            eta_s[replace_idx] = new_eta_s
        end
        
        # Update η population (variable to clause messages)
        for i in 1:N
            p = rand(rng, Poisson(k * alpha / 2))
            q = rand(rng, Poisson(k * alpha / 2))
            total = p + q
            
            if total == 0
                new_eta = 0.5
            else
                selected = rand(rng, 1:N, total)
                pos_neg = eta_s[selected]
                pos = @view pos_neg[1:p]
                neg = @view pos_neg[p+1:end]
                
                # Stabilize computation using logarithms
                log_prod1 = sum(log, pos) + sum(x -> log(1 - x), neg)
                log_prod2 = sum(x -> log(1 - x), pos) + sum(log, neg)
                max_log = max(log_prod1, log_prod2)
                log_total = max_log + log(exp(log_prod1 - max_log) + exp(log_prod2 - max_log))
                new_eta = exp(log_prod1 - log_total)
            end
            
            replace_idx = rand(rng, 1:N)
            eta[replace_idx] = new_eta
        end
    end
    
    # Calculate entropy components
    term1 = zeros(num_samples)
    term2 = zeros(num_samples)
    term3 = zeros(num_samples)
    
    # Parallel computation using threads
    Threads.@threads for i in 1:num_samples
        rng = Random.MersenneTwister()  # Thread-local RNG
        
        # Term 1 calculation
        indices = rand(rng, 1:N, k)
        product = prod(eta[indices])
        term1[i] = log(1 - product + 1e-10)
        
        # Term 2 calculation
        p = rand(rng, Poisson(k * alpha / 2))
        q = rand(rng, Poisson(k * alpha / 2))
        total = p + q
        
        if total == 0
            term2[i] = log(2.0)
        else
            selected = rand(rng, 1:N, total)
            pos_neg = eta_s[selected]
            pos = @view pos_neg[1:p]
            neg = @view pos_neg[p+1:end]
            
            log_prod1 = sum(log, pos) + sum(x -> log(1 - x), neg)
            log_prod2 = sum(x -> log(1 - x), pos) + sum(log, neg)
            max_log = max(log_prod1, log_prod2)
            term2[i] = max_log + log(exp(log_prod1 - max_log) + exp(log_prod2 - max_log))
        end
        
        # Term 3 calculation
        idx_eta = rand(rng, 1:N)
        idx_eta_s = rand(rng, 1:N)
        term3[i] = log(eta_s[idx_eta_s] * eta[idx_eta] + 
                      (1 - eta_s[idx_eta_s]) * (1 - eta[idx_eta]) + 1e-10)
    end
    
    # Final entropy calculation
    s = alpha * mean(term1) + mean(term2) - (k * alpha / 2) * mean(term3)
    return s
end

# Example usage
let alpha = 4.0  # Clause density
    entropy = compute_entropy(alpha)
    println("Entropy density: ", entropy)
end