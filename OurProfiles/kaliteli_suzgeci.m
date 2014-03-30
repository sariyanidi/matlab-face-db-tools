function kaliteli_suzgeci
% Profil yuzler icin isaretledigimiz yuzlerden sadece buyuk olanları (48x48
% piksel+ustu) tutacak olan bir suzgec dosyasi
clc

girdi_dizini = '../yeni_koordinatlar';
cikti_dizini = '../suzulen_imajlar';

% yuzlerin orjinal resimlere gore isaretlenmis karelerini iceren dosyalar
koord_dosyalari = cell(5,1);
koord_dosyalari{1} = 'yeni_koordinatlar_BIRKAN.txt';
koord_dosyalari{2} = 'yeni_koordinatlar_cihan.txt';
koord_dosyalari{3} = 'yeni_koordinatlar_bora.txt';
koord_dosyalari{4} = 'yeni_koordinatlar_VANGEL.txt';
koord_dosyalari{5} = 'yeni_koordinatlar_VOLKAN.txt';

for i=1:length(koord_dosyalari)
    dosya_sahibi = koord_dosyalari{i}(19:end-4);
    taze_dosya = fopen([cikti_dizini '/' dosya_sahibi '.txt'], 'w');
    
    fprintf('%s efendinin kesimleri irdeleniyor \n', dosya_sahibi);
    %fprintf('========================================= \n');
    
    girdi_yolu = [girdi_dizini '/' koord_dosyalari{i}];
    fid = fopen(girdi_yolu, 'r');
    
    if -1 == fid
        disp('Text dosyasi acilamadi! devam et...');
        return;
    end
    
    j=-1;
    gecerli_yuz_sayisi = 0;
    dosya_satiri = fgetl(fid);
    
    while dosya_satiri ~= -1
        j = j+1;
        %fprintf('%d... ', j);
        temp = textscan(dosya_satiri,'%s\t%d\t%d\t%d\t%d\t%d',1);
        ham_dosya_ismi = cell2mat(temp{1});
        temp = cell2mat(temp(2:6));
        poz = temp(5);
        dikd = temp(1:4);
        dosya_satiri = fgetl(fid);
        
        if poz == 3
            j=j-1;
            continue;
        end
        
        if (temp(3) < 100)
            %fprintf('yuz cok kucuk! (%d piksel)\n', temp(3));
            continue;
        end
        
        
        fprintf(taze_dosya, 'image%d %s %d %d %d %d\n', j, ham_dosya_ismi, dikd(1), dikd(2), dikd(3), dikd(4));
        %fprintf('yuz eklendi ;) (%d piksel boyunda)\n', temp(3));
        
        gecerli_yuz_sayisi = gecerli_yuz_sayisi+1;
    end
    
    if gecerli_yuz_sayisi < 300
        fprintf('%s! Topu topu %d adam akilli yuz cikmis senden... Yaziklar olsun, bu projede yerin yok...\n', dosya_sahibi, gecerli_yuz_sayisi);
    elseif gecerli_yuz_sayisi > 700
        fprintf('%s abi, masallah %d tane on numara yuz secmissin, ben bir yazilim olarak en kaliteli kesimlerin seninkiler olduguna karar verdim.\n', dosya_sahibi, gecerli_yuz_sayisi);
        fprintf('Basarilarinin devamini diliyorum abi.\n');
    else
        fprintf('Afferin lan %s, az bucuk duzgun biseyler kesmissin bari sen, %d tane duzgun yuz var.\n', dosya_sahibi, gecerli_yuz_sayisi);
    end
    
    fclose(taze_dosya);
    fprintf('--------------------------------------------------------------------\n\n');
end


end