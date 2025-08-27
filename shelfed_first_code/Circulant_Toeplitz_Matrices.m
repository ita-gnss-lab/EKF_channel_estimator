function Circulant_Matrices = Circulant_Toeplitz_Matrices(Size, Quantity)
    
    Circulant_Matrices = [];

    for Potency = 0 : Quantity
        Circulant_Matrices(:, :, Potency+1) = [zeros(Potency, Size - Potency) eye(Potency);
                                               eye(Size - Potency) zeros(Size - Potency, Potency)];
    end

end
