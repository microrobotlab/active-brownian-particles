# Functional forms of most used forces in ABP main

lennard_jones(x, σ, ϵ) = 24*ϵ*(((2*σ^(12))/(x^(13)))- (σ^(6)/(x^(7))))
shifted_lennard_jones(x, σ, ϵ, shift) = 24*ϵ*(((2*σ^(12))/((x-shift)^(13)))- (σ^(6)/((x-shift)^(7))))

function lj_nondiv(x, σ, ϵ)
    if x > σ
        return lennard_jones(x, σ, ϵ)
    else
        return lennard_jones(σ, σ, ϵ)
    end    
end

function weeks_chandler_andersen(x, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6)
    if x > rmin
        return 0
    else
        return 24*ϵ*(((2*σ^(12))/(x^(13)))- (σ^(6)/(x^(7))))
    end
end

function purely_attractive_lennard_jones(x, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 24*ϵ*(((2*σ^(12))/((x + σ/2)^(13)))- (σ^(6)/((x + σ/2)^(7))))
    else
        return 0
    end
end

function rlj_boundary(x::Float64, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 0

    elseif x < σ
        return 390144*ϵ/σ
    else
        return 24*ϵ*(((2*σ^(12))/((x-σ/2)^(13)))- (σ^(6)/((x-σ/2)^(7)))) 
    end
end

function rlj_boundary(x::Float64, σ::Float64, ϵ::Float64)
    rmin = σ*2^(1/6) + σ/2
    if x > rmin
        return 0

    elseif x < σ
        return 390144*ϵ/σ
    else
        return 24*ϵ*(((2*σ^(12))/((x-σ/2)^(13)))- (σ^(6)/((x-σ/2)^(7)))) 
    end
end

function contact_lj(x, σ::Float64, ϵ::Float64)
    shift = 2^(1/6)-1
    return 24*ϵ*(((2*σ^(12))/((x+shift)^(13)))- (σ^(6)/((x+shift)^(7))))
end

coulomb(x,k) = k / x^2

excluded_volume(x::Float64,R::Float64, ϵ::Float64) = 12ϵ*((2R)^12)x^(-13)

spring(x,k, x0) = -k*(x-x0)