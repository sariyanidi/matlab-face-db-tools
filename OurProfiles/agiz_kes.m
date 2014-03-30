function agiz_kes
% Kendi topladigimiz profil yuzleri ve isaretledigimiz koordinatlari
% kullanarak otomatik agiz kesmece yap

ana_girdi_dizini = '/media/DEPO/Server_Yedek_20110512/veritabani/Saptama/acili_yuzler/HizaliKoordinatlar';
suzulmus_dizini = '/media/DEPO/Server_Yedek_20110512/veritabani/Saptama/acili_yuzler/suzulen_imajlar';
ana_cikti_dizini = '/media/DEPO/Server_Yedek_20110512/veritabani/Saptama/acili_yuzler/agiz_kesilmis';

kisiler = {'birkan', 'bora', 'cihan', 'vangel','volkan'};
  
for k=1:length(kisiler)
    kisi_adi = kisiler{k};
    fprintf('%s adli sahsiyetin dosyalarindan agiz kesiyoruz\n', kisi_adi);
    fprintf('-------------------------------------------------------\n');
    
    suzulmus_dosya = fopen([suzulmus_dizini '/' kisi_adi '.txt'], 'r');
    
    imaj_dizini = [ana_girdi_dizini '/imageFiles_' kisi_adi];
    koor_dizini = [ana_girdi_dizini '/coordinateFiles_' kisi_adi];
    cikti_dizini = [ana_cikti_dizini '/' kisi_adi];
    
    if ~isdir(cikti_dizini)
        mkdir(cikti_dizini)
    end
    
    suzulmus_satir = fgetl(suzulmus_dosya);
    toplam_kesilen = 0;
    
    while suzulmus_satir ~= -1
        temp = textscan(suzulmus_satir,'%s\t%s\t%d\t%d\t%d\t%d',1);
        suzulmus_satir = fgetl(suzulmus_dosya);
        
        hedef_imaj = dir([imaj_dizini '/' cell2mat(temp{1}) '.png']);
        hedef_koor = dir([koor_dizini '/' cell2mat(temp{1}) '__*']);
        
        % dosya bulunamiyorsa vs.. devam et.
        try
            idx = strfind(hedef_koor.name, '__');
            poz = str2double(hedef_koor.name(idx+5));
            
        im = imread([imaj_dizini '/' hedef_imaj.name]);
        nokta = load([koor_dizini '/' hedef_koor.name]);
        catch e
            fprintf('Dikkat, hata: %s (dosya adi: %s)\n', e.message, cell2mat(temp{1}));
            continue;
        end
        
        % poz33'lu dosya isimlerini boyle ayikla
        if 3 == str2double(hedef_koor.name(idx+6))
            continue;
        end
        
        
        width = size(im,2);
        height = size(im,1);
        if poz == 4 || poz == 5
            im = flipdim(im,2);
            nokta =  nokta_fliple(nokta, width);
            poz = 6-poz;
        end
        
        % Profil ise
        if poz == 1
            aci = atan((nokta(2,2)-nokta(4,2)) / (nokta(4,1)-nokta(2,1)));
            nokta = nokta_dondur(nokta', aci, width, height);
            
            % noktalar dondu, burun ve ceneyi tekrar tanimla
            burun = nokta(2,:);
            cene = nokta(3,:);
            % Ara aci ise
        elseif poz == 2 
            burun = nokta(3,:);
            temp = [nokta(1,1), nokta(2,1)];
            n1 = nokta(temp == min(temp), :);
            n2 = nokta(temp == max(temp), :);
            
            aci = atan((n1(2)-n2(2)) / (n2(1)-n1(1)));
            
            % iki goz birbirine cok yakin olunca d.ici donme degerleri cok
            % hatali geliyor, bu durumlarda duzlem ici dondurme burun ile
            % kulak arasindaki acidan aliniyor
            if abs(n1(1)-n2(1)) < 18
                kulak = nokta(5,:);
                aci = -atan((burun(2)-kulak(2))/(burun(1)-kulak(1)));
            end
            nokta = nokta_dondur(nokta', aci, width, height);
            
            % noktalar dondu, burun ve ceneyi tekrar tanimla
            burun = nokta(3,:);
            cene = nokta(4,:);
        end
        
        % yuzu bi dondur hele
        aci_deg = (180 * -aci) / pi;
        im = imrotate(im, aci_deg, 'bicubic', 'crop');
        
        % agiz karesini belirle
        agiz_g = cene(2)-burun(2);

        if 1 == poz
            slide_offset = 0.15;
        elseif 2 == poz
            slide_offset = 0.1;
        end
        
        % burun ile cenenin tam ortasi - agiz karesinin x koordinati buraya
        % gore belirenecek
        orta_nokta = [0.5*(burun(1)+cene(1)); 
                      0.5*(burun(2)+cene(2))];
        
        agiz_karesi = [orta_nokta(1)-0.5*agiz_g;
                       burun(2)+0.0*agiz_g;
                       1.0*agiz_g;
                       1.0*agiz_g];
                   
        agiz_karesi(1) = agiz_karesi(1)+agiz_karesi(1)*slide_offset;

        % agizi kes
        agiz = imcrop(im, agiz_karesi);
        if 3 == length(size(agiz))
            agiz = rgb2gray(agiz);
        end
        
        % cok kucuk agizlari ele
        if 20 > size(agiz,1); continue; end
        imwrite(agiz, [cikti_dizini '/' hedef_imaj.name]);
        
        toplam_kesilen = toplam_kesilen+1;
        if mod(toplam_kesilen,100) == 0
            fprintf('%d nolu agiz da kesildi\n', toplam_kesilen);
        end
    end
    
    fprintf('----------------------------------------------------\n\n');
    fclose(suzulmus_dosya);
end

end