function flipli = nokta_fliple(nokta, iw)
    nokta(:,1) = iw-nokta(:,1);
    flipli = nokta;
end