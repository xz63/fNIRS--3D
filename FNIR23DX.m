classdef FNIR23DX <hgsetget
    %  this program is mde by Dr. Xian Zhang fom brain function lab, Yale School of Medicine.
    properties
        polhemus % used for saving the result, reading the original data ,,,
        img
        %imgSkull
        xi
        yi
        zi
        x
        y
        z
        V
        w
        M
    end % end of propertieso
    methods
        function x = FNIR23DX ()
            load('MNITemplate');
            d=9;
            m=zeros(d,d,d);
            m=m+1;
            img1=imerode(img,m);
            img1=img-img1;
            x.img=img1;
            %load('MNITemplate_skull');
            %x.imgSkull=img;
            sz=size(img);
            load MNITemplate_Header
            x.V=V;
        end
        function x=setXYZ(x,xyz);
            cor = mni2cor(xyz, x.V.mat);
            x.x=cor(:,1);
            x.y=cor(:,2);
            x.z=cor(:,3);
            [x.xi,x.yi,x.zi] = meshgrid([(-4+min(cor(:,1))):(4+max(cor(:,1)))], [(-4+min(cor(:,2))):(+4+max(cor(:,2)))], [(-4+min(cor(:,3))):(4+max(cor(:,3)))]);
        end
        function x=setPolhemus(x,polhemus);
            x.polhemus =polhemus;
            nCh=length(x.polhemus.chname);
            xyz=x.polhemus.NFRI_result.OtherC(x.polhemus.nTR+[1:nCh],:);
            x.setXYZ(xyz);
        end
        
        
        function M=getValue(x,beta,fn);
            x.V.fname=[fn '.nii'];
            %            if exist(x.V.fname); return; end
            %            fid=fopen(x.V.fname,'w');
            %            fprintf(fid,'s');
            %            fclose(fid);
            for i=1:size(beta,2);
                v = griddatan([x.x(:) x.y(:) x.z(:)], beta(:,i), [x.xi(:)  x.yi(:)  x.zi(:)] );
                M=x.img*NaN;
                for ii=1:length(x.xi(:))
                    if (~isnan(v(ii))) &  (x.img(x.xi(ii),x.yi(ii),x.zi(ii))>0)
                        M(x.xi(ii),x.yi(ii),x.zi(ii),i)= v(ii);
                    end
                end
            end
        end
        function x=save(x,fn,M);
            x.V.fname=[fn '.nii'];
	    if exist(x.V.fname); return;end
            for i=1:size(M,4);
                spm_write_vol(x.V,M(:,:,:,i));
            end
        end
        function renderChannel_zscore(x,fn,z)
            p1(:,1)=x.xi(:);
            p1(:,2)=x.yi(:);
            p1(:,3)=x.zi(:);
            img=x.img*0;
            for ch=1:length(z)
                p2(:,1)=x.x(ch);
                p2(:,2)=x.y(ch);
                p2(:,3)=x.z(ch);
                pd=pdist2(p1,p2);
                pd1=min(pd,[],2);
                for i=1:length(x.xi(:))
                    %img(x.xi(i),x.yi(i),x.zi(i))=z(ch)*exp( -pd1(i)/2);
                    %if pd1(i)<3; img(x.xi(i),x.yi(i),x.zi(i))=z(ch);end
                    if pd1(i)<15; img(x.xi(i),x.yi(i),x.zi(i),ch)=z(ch)*exp(-(pd1(i)).^2/30);end
                end
            end
            img=sum(img,4);
            x.V.fname=[fn '.nii'];
            spm_write_vol(x.V,img);
        end
        
        function renderChannel(x,fn)
            p1(:,1)=x.xi(:);
            p1(:,2)=x.yi(:);
            p1(:,3)=x.zi(:);
            p2(:,1)=x.x;
            p2(:,2)=x.y;
            p2(:,3)=x.z;
            
            pd=pdist2(p1,p2);
            pd1=min(pd,[],2);
            img=x.img*NaN;
            for i=1:length(x.xi(:))
                img(x.xi(i),x.yi(i),x.zi(i))=10*exp( -pd1(i)/2);
            end
            x.V.fname=[fn '.nii'];
            spm_write_vol(x.V,img);
        end
        function demo() % you need copy and paste and make a script
            rmsv=squeeze(mean(rms(cleandata(:,1,:,1,:)),5));
            f=FNIR23DX();
            %            f.setPolhemus(SubjectGroupChannel{1}.polhemus);
            f.setXYZ(xyz3d);
            M=f.getValue(rmsv,'result');
            f.save('result',M);
        end
        function M=renderChannelValue(x,value, r,fn)
            v = griddatan([x.x(:) x.y(:) x.z(:)], value, [x.xi(:)  x.yi(:)  x.zi(:)] , 'nearest');
            M=x.img*NaN;
            for ii=1:length(x.xi(:))
                if (~isnan(v(ii))) &  (x.img(x.xi(ii),x.yi(ii),x.zi(ii))>0)
                    M(x.xi(ii),x.yi(ii),x.zi(ii))= v(ii);
                end
            end
            p1(:,1)=x.xi(:);
            p1(:,2)=x.yi(:);
            p1(:,3)=x.zi(:);
            p2(:,1)=x.x;
            p2(:,2)=x.y;
            p2(:,3)=x.z;
            
            pd=pdist2(p1,p2);
            pd1=min(pd,[],2);
            img=x.img*NaN;
            for i=1:length(x.xi(:))
                img(x.xi(i),x.yi(i),x.zi(i))=(r>pd1(i))*exp( -pd1(i)/2);
            end
            M=img.*M;
            x.V.fname=[fn '.nii'];
            spm_write_vol(x.V,M);
        end
    end
end
