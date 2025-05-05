function [B_prime, Beta_g] = vsfp3_gsr(input_data,hem_data)

    B = input_data';
    g = hem_data';
    g_plus = (1 / (g' * g)) * g';
    Beta_g = g_plus * B;
    B_prime = B - g * Beta_g;      
end