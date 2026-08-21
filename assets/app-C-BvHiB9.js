const __vite__mapDeps=(i,m=__vite__mapDeps,d=(m.f||(m.f=["assets/SearchModal-BnGoJJQr.js","assets/SearchModal-Cu0-bXnC.css"])))=>i.map(i=>d[i]);
(function(){const n=document.createElement("link").relList;if(n&&n.supports&&n.supports("modulepreload"))return;for(const t of document.querySelectorAll('link[rel="modulepreload"]'))e(t);new MutationObserver(t=>{for(const l of t)if(l.type==="childList")for(const o of l.addedNodes)o.tagName==="LINK"&&o.rel==="modulepreload"&&e(o)}).observe(document,{childList:!0,subtree:!0});function a(t){const l={};return t.integrity&&(l.integrity=t.integrity),t.referrerPolicy&&(l.referrerPolicy=t.referrerPolicy),t.crossOrigin==="use-credentials"?l.credentials="include":t.crossOrigin==="anonymous"?l.credentials="omit":l.credentials="same-origin",l}function e(t){if(t.ep)return;t.ep=!0;const l=a(t);fetch(t.href,l)}})();const Iu="modulepreload",Nu=function(s){return"/"+s},Ko={},Cs=function(n,a,e){let t=Promise.resolve();if(a&&a.length>0){let c=function(r){return Promise.all(r.map(i=>Promise.resolve(i).then(u=>({status:"fulfilled",value:u}),u=>({status:"rejected",reason:u}))))};document.getElementsByTagName("link");const o=document.querySelector("meta[property=csp-nonce]"),p=o?.nonce||o?.getAttribute("nonce");t=c(a.map(r=>{if(r=Nu(r),r in Ko)return;Ko[r]=!0;const i=r.endsWith(".css"),u=i?'[rel="stylesheet"]':"";if(document.querySelector(`link[href="${r}"]${u}`))return;const h=document.createElement("link");if(h.rel=i?"stylesheet":Iu,i||(h.as="script"),h.crossOrigin="",h.href=r,p&&h.setAttribute("nonce",p),document.head.appendChild(h),i)return new Promise((d,f)=>{h.addEventListener("load",d),h.addEventListener("error",()=>f(new Error(`Unable to preload CSS for ${r}`)))})}))}function l(o){const p=new Event("vite:preloadError",{cancelable:!0});if(p.payload=o,window.dispatchEvent(p),!p.defaultPrevented)throw o}return t.then(o=>{for(const p of o||[])p.status==="rejected"&&l(p.reason);return n().catch(l)})};/**
* @vue/shared v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function to(s){const n=Object.create(null);for(const a of s.split(","))n[a]=1;return a=>a in n}const As={},Wa=[],qn=()=>{},Wc=()=>!1,Rt=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&(s.charCodeAt(2)>122||s.charCodeAt(2)<97),lo=s=>s.startsWith("onUpdate:"),Ys=Object.assign,oo=(s,n)=>{const a=s.indexOf(n);a>-1&&s.splice(a,1)},Fu=Object.prototype.hasOwnProperty,Ss=(s,n)=>Fu.call(s,n),rs=Array.isArray,Ka=s=>Lt(s)==="[object Map]",Kc=s=>Lt(s)==="[object Set]",ds=s=>typeof s=="function",Us=s=>typeof s=="string",la=s=>typeof s=="symbol",Ms=s=>s!==null&&typeof s=="object",Xc=s=>(Ms(s)||ds(s))&&ds(s.then)&&ds(s.catch),Yc=Object.prototype.toString,Lt=s=>Yc.call(s),$u=s=>Lt(s).slice(8,-1),Qc=s=>Lt(s)==="[object Object]",po=s=>Us(s)&&s!=="NaN"&&s[0]!=="-"&&""+parseInt(s,10)===s,de=to(",key,ref,ref_for,ref_key,onVnodeBeforeMount,onVnodeMounted,onVnodeBeforeUpdate,onVnodeUpdated,onVnodeBeforeUnmount,onVnodeUnmounted"),Mt=s=>{const n=Object.create(null);return(a=>n[a]||(n[a]=s(a)))},Bu=/-\w/g,Tn=Mt(s=>s.replace(Bu,n=>n.slice(1).toUpperCase())),qu=/\B([A-Z])/g,Ca=Mt(s=>s.replace(qu,"-$1").toLowerCase()),Dt=Mt(s=>s.charAt(0).toUpperCase()+s.slice(1)),Qt=Mt(s=>s?`on${Dt(s)}`:""),ba=(s,n)=>!Object.is(s,n),ot=(s,...n)=>{for(let a=0;a<s.length;a++)s[a](...n)},Jc=(s,n,a,e=!1)=>{Object.defineProperty(s,n,{configurable:!0,enumerable:!1,writable:e,value:a})},co=s=>{const n=parseFloat(s);return isNaN(n)?s:n},zu=s=>{const n=Us(s)?Number(s):NaN;return isNaN(n)?s:n};let Xo;const Ot=()=>Xo||(Xo=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:{});function oa(s){if(rs(s)){const n={};for(let a=0;a<s.length;a++){const e=s[a],t=Us(e)?Gu(e):oa(e);if(t)for(const l in t)n[l]=t[l]}return n}else if(Us(s)||Ms(s))return s}const Uu=/;(?![^(]*\))/g,Hu=/:([^]+)/,Vu=/\/\*[^]*?\*\//g;function Gu(s){const n={};return s.replace(Vu,"").split(Uu).forEach(a=>{if(a){const e=a.split(Hu);e.length>1&&(n[e[0].trim()]=e[1].trim())}}),n}function Fs(s){let n="";if(Us(s))n=s;else if(rs(s))for(let a=0;a<s.length;a++){const e=Fs(s[a]);e&&(n+=e+" ")}else if(Ms(s))for(const a in s)s[a]&&(n+=a+" ");return n.trim()}const Wu="itemscope,allowfullscreen,formnovalidate,ismap,nomodule,novalidate,readonly",Ku=to(Wu);function Zc(s){return!!s||s===""}const sr=s=>!!(s&&s.__v_isRef===!0),V=s=>Us(s)?s:s==null?"":rs(s)||Ms(s)&&(s.toString===Yc||!ds(s.toString))?sr(s)?V(s.value):JSON.stringify(s,nr,2):String(s),nr=(s,n)=>sr(n)?nr(s,n.value):Ka(n)?{[`Map(${n.size})`]:[...n.entries()].reduce((a,[e,t],l)=>(a[Jt(e,l)+" =>"]=t,a),{})}:Kc(n)?{[`Set(${n.size})`]:[...n.values()].map(a=>Jt(a))}:la(n)?Jt(n):Ms(n)&&!rs(n)&&!Qc(n)?String(n):n,Jt=(s,n="")=>{var a;return la(s)?`Symbol(${(a=s.description)!=null?a:n})`:s};/**
* @vue/reactivity v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let ln;class ar{constructor(n=!1){this.detached=n,this._active=!0,this._on=0,this.effects=[],this.cleanups=[],this._isPaused=!1,this.parent=ln,!n&&ln&&(this.index=(ln.scopes||(ln.scopes=[])).push(this)-1)}get active(){return this._active}pause(){if(this._active){this._isPaused=!0;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].pause();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].pause()}}resume(){if(this._active&&this._isPaused){this._isPaused=!1;let n,a;if(this.scopes)for(n=0,a=this.scopes.length;n<a;n++)this.scopes[n].resume();for(n=0,a=this.effects.length;n<a;n++)this.effects[n].resume()}}run(n){if(this._active){const a=ln;try{return ln=this,n()}finally{ln=a}}}on(){++this._on===1&&(this.prevScope=ln,ln=this)}off(){this._on>0&&--this._on===0&&(ln=this.prevScope,this.prevScope=void 0)}stop(n){if(this._active){this._active=!1;let a,e;for(a=0,e=this.effects.length;a<e;a++)this.effects[a].stop();for(this.effects.length=0,a=0,e=this.cleanups.length;a<e;a++)this.cleanups[a]();if(this.cleanups.length=0,this.scopes){for(a=0,e=this.scopes.length;a<e;a++)this.scopes[a].stop(!0);this.scopes.length=0}if(!this.detached&&this.parent&&!n){const t=this.parent.scopes.pop();t&&t!==this&&(this.parent.scopes[this.index]=t,t.index=this.index)}this.parent=void 0}}}function ro(s){return new ar(s)}function er(){return ln}function io(s,n=!1){ln&&ln.cleanups.push(s)}let Ls;const Zt=new WeakSet;class tr{constructor(n){this.fn=n,this.deps=void 0,this.depsTail=void 0,this.flags=5,this.next=void 0,this.cleanup=void 0,this.scheduler=void 0,ln&&ln.active&&ln.effects.push(this)}pause(){this.flags|=64}resume(){this.flags&64&&(this.flags&=-65,Zt.has(this)&&(Zt.delete(this),this.trigger()))}notify(){this.flags&2&&!(this.flags&32)||this.flags&8||or(this)}run(){if(!(this.flags&1))return this.fn();this.flags|=2,Yo(this),pr(this);const n=Ls,a=Rn;Ls=this,Rn=!0;try{return this.fn()}finally{cr(this),Ls=n,Rn=a,this.flags&=-3}}stop(){if(this.flags&1){for(let n=this.deps;n;n=n.nextDep)mo(n);this.deps=this.depsTail=void 0,Yo(this),this.onStop&&this.onStop(),this.flags&=-2}}trigger(){this.flags&64?Zt.add(this):this.scheduler?this.scheduler():this.runIfDirty()}runIfDirty(){wl(this)&&this.run()}get dirty(){return wl(this)}}let lr=0,me,je;function or(s,n=!1){if(s.flags|=8,n){s.next=je,je=s;return}s.next=me,me=s}function uo(){lr++}function ho(){if(--lr>0)return;if(je){let n=je;for(je=void 0;n;){const a=n.next;n.next=void 0,n.flags&=-9,n=a}}let s;for(;me;){let n=me;for(me=void 0;n;){const a=n.next;if(n.next=void 0,n.flags&=-9,n.flags&1)try{n.trigger()}catch(e){s||(s=e)}n=a}}if(s)throw s}function pr(s){for(let n=s.deps;n;n=n.nextDep)n.version=-1,n.prevActiveLink=n.dep.activeLink,n.dep.activeLink=n}function cr(s){let n,a=s.depsTail,e=a;for(;e;){const t=e.prevDep;e.version===-1?(e===a&&(a=t),mo(e),Xu(e)):n=e,e.dep.activeLink=e.prevActiveLink,e.prevActiveLink=void 0,e=t}s.deps=n,s.depsTail=a}function wl(s){for(let n=s.deps;n;n=n.nextDep)if(n.dep.version!==n.version||n.dep.computed&&(rr(n.dep.computed)||n.dep.version!==n.version))return!0;return!!s._dirty}function rr(s){if(s.flags&4&&!(s.flags&16)||(s.flags&=-17,s.globalVersion===ke)||(s.globalVersion=ke,!s.isSSR&&s.flags&128&&(!s.deps&&!s._dirty||!wl(s))))return;s.flags|=2;const n=s.dep,a=Ls,e=Rn;Ls=s,Rn=!0;try{pr(s);const t=s.fn(s._value);(n.version===0||ba(t,s._value))&&(s.flags|=128,s._value=t,n.version++)}catch(t){throw n.version++,t}finally{Ls=a,Rn=e,cr(s),s.flags&=-3}}function mo(s,n=!1){const{dep:a,prevSub:e,nextSub:t}=s;if(e&&(e.nextSub=t,s.prevSub=void 0),t&&(t.prevSub=e,s.nextSub=void 0),a.subs===s&&(a.subs=e,!e&&a.computed)){a.computed.flags&=-5;for(let l=a.computed.deps;l;l=l.nextDep)mo(l,!0)}!n&&!--a.sc&&a.map&&a.map.delete(a.key)}function Xu(s){const{prevDep:n,nextDep:a}=s;n&&(n.nextDep=a,s.prevDep=void 0),a&&(a.prevDep=n,s.nextDep=void 0)}let Rn=!0;const ir=[];function sa(){ir.push(Rn),Rn=!1}function na(){const s=ir.pop();Rn=s===void 0?!0:s}function Yo(s){const{cleanup:n}=s;if(s.cleanup=void 0,n){const a=Ls;Ls=void 0;try{n()}finally{Ls=a}}}let ke=0;class Yu{constructor(n,a){this.sub=n,this.dep=a,this.version=a.version,this.nextDep=this.prevDep=this.nextSub=this.prevSub=this.prevActiveLink=void 0}}class jo{constructor(n){this.computed=n,this.version=0,this.activeLink=void 0,this.subs=void 0,this.map=void 0,this.key=void 0,this.sc=0,this.__v_skip=!0}track(n){if(!Ls||!Rn||Ls===this.computed)return;let a=this.activeLink;if(a===void 0||a.sub!==Ls)a=this.activeLink=new Yu(Ls,this),Ls.deps?(a.prevDep=Ls.depsTail,Ls.depsTail.nextDep=a,Ls.depsTail=a):Ls.deps=Ls.depsTail=a,ur(a);else if(a.version===-1&&(a.version=this.version,a.nextDep)){const e=a.nextDep;e.prevDep=a.prevDep,a.prevDep&&(a.prevDep.nextDep=e),a.prevDep=Ls.depsTail,a.nextDep=void 0,Ls.depsTail.nextDep=a,Ls.depsTail=a,Ls.deps===a&&(Ls.deps=e)}return a}trigger(n){this.version++,ke++,this.notify(n)}notify(n){uo();try{for(let a=this.subs;a;a=a.prevSub)a.sub.notify()&&a.sub.dep.notify()}finally{ho()}}}function ur(s){if(s.dep.sc++,s.sub.flags&4){const n=s.dep.computed;if(n&&!s.dep.subs){n.flags|=20;for(let e=n.deps;e;e=e.nextDep)ur(e)}const a=s.dep.subs;a!==s&&(s.prevSub=a,a&&(a.nextSub=s)),s.dep.subs=s}}const mt=new WeakMap,Ma=Symbol(""),Cl=Symbol(""),xe=Symbol("");function pn(s,n,a){if(Rn&&Ls){let e=mt.get(s);e||mt.set(s,e=new Map);let t=e.get(a);t||(e.set(a,t=new jo),t.map=e,t.key=a),t.track()}}function Yn(s,n,a,e,t,l){const o=mt.get(s);if(!o){ke++;return}const p=c=>{c&&c.trigger()};if(uo(),n==="clear")o.forEach(p);else{const c=rs(s),r=c&&po(a);if(c&&a==="length"){const i=Number(e);o.forEach((u,h)=>{(h==="length"||h===xe||!la(h)&&h>=i)&&p(u)})}else switch((a!==void 0||o.has(void 0))&&p(o.get(a)),r&&p(o.get(xe)),n){case"add":c?r&&p(o.get("length")):(p(o.get(Ma)),Ka(s)&&p(o.get(Cl)));break;case"delete":c||(p(o.get(Ma)),Ka(s)&&p(o.get(Cl)));break;case"set":Ka(s)&&p(o.get(Ma));break}}ho()}function Qu(s,n){const a=mt.get(s);return a&&a.get(n)}function $a(s){const n=ws(s);return n===s?n:(pn(n,"iterate",xe),Sn(s)?n:n.map(Zs))}function It(s){return pn(s=ws(s),"iterate",xe),s}const Ju={__proto__:null,[Symbol.iterator](){return sl(this,Symbol.iterator,Zs)},concat(...s){return $a(this).concat(...s.map(n=>rs(n)?$a(n):n))},entries(){return sl(this,"entries",s=>(s[1]=Zs(s[1]),s))},every(s,n){return Hn(this,"every",s,n,void 0,arguments)},filter(s,n){return Hn(this,"filter",s,n,a=>a.map(Zs),arguments)},find(s,n){return Hn(this,"find",s,n,Zs,arguments)},findIndex(s,n){return Hn(this,"findIndex",s,n,void 0,arguments)},findLast(s,n){return Hn(this,"findLast",s,n,Zs,arguments)},findLastIndex(s,n){return Hn(this,"findLastIndex",s,n,void 0,arguments)},forEach(s,n){return Hn(this,"forEach",s,n,void 0,arguments)},includes(...s){return nl(this,"includes",s)},indexOf(...s){return nl(this,"indexOf",s)},join(s){return $a(this).join(s)},lastIndexOf(...s){return nl(this,"lastIndexOf",s)},map(s,n){return Hn(this,"map",s,n,void 0,arguments)},pop(){return oe(this,"pop")},push(...s){return oe(this,"push",s)},reduce(s,...n){return Qo(this,"reduce",s,n)},reduceRight(s,...n){return Qo(this,"reduceRight",s,n)},shift(){return oe(this,"shift")},some(s,n){return Hn(this,"some",s,n,void 0,arguments)},splice(...s){return oe(this,"splice",s)},toReversed(){return $a(this).toReversed()},toSorted(s){return $a(this).toSorted(s)},toSpliced(...s){return $a(this).toSpliced(...s)},unshift(...s){return oe(this,"unshift",s)},values(){return sl(this,"values",Zs)}};function sl(s,n,a){const e=It(s),t=e[n]();return e!==s&&!Sn(s)&&(t._next=t.next,t.next=()=>{const l=t._next();return l.done||(l.value=a(l.value)),l}),t}const Zu=Array.prototype;function Hn(s,n,a,e,t,l){const o=It(s),p=o!==s&&!Sn(s),c=o[n];if(c!==Zu[n]){const u=c.apply(s,l);return p?Zs(u):u}let r=a;o!==s&&(p?r=function(u,h){return a.call(this,Zs(u),h,s)}:a.length>2&&(r=function(u,h){return a.call(this,u,h,s)}));const i=c.call(o,r,e);return p&&t?t(i):i}function Qo(s,n,a,e){const t=It(s);let l=a;return t!==s&&(Sn(s)?a.length>3&&(l=function(o,p,c){return a.call(this,o,p,c,s)}):l=function(o,p,c){return a.call(this,o,Zs(p),c,s)}),t[n](l,...e)}function nl(s,n,a){const e=ws(s);pn(e,"iterate",xe);const t=e[n](...a);return(t===-1||t===!1)&&bo(a[0])?(a[0]=ws(a[0]),e[n](...a)):t}function oe(s,n,a=[]){sa(),uo();const e=ws(s)[n].apply(s,a);return ho(),na(),e}const sh=to("__proto__,__v_isRef,__isVue"),hr=new Set(Object.getOwnPropertyNames(Symbol).filter(s=>s!=="arguments"&&s!=="caller").map(s=>Symbol[s]).filter(la));function nh(s){la(s)||(s=String(s));const n=ws(this);return pn(n,"has",s),n.hasOwnProperty(s)}class dr{constructor(n=!1,a=!1){this._isReadonly=n,this._isShallow=a}get(n,a,e){if(a==="__v_skip")return n.__v_skip;const t=this._isReadonly,l=this._isShallow;if(a==="__v_isReactive")return!t;if(a==="__v_isReadonly")return t;if(a==="__v_isShallow")return l;if(a==="__v_raw")return e===(t?l?uh:gr:l?fr:jr).get(n)||Object.getPrototypeOf(n)===Object.getPrototypeOf(e)?n:void 0;const o=rs(n);if(!t){let c;if(o&&(c=Ju[a]))return c;if(a==="hasOwnProperty")return nh}const p=Reflect.get(n,a,zs(n)?n:e);if((la(a)?hr.has(a):sh(a))||(t||pn(n,"get",a),l))return p;if(zs(p)){const c=o&&po(a)?p:p.value;return t&&Ms(c)?xl(c):c}return Ms(p)?t?xl(p):$e(p):p}}class mr extends dr{constructor(n=!1){super(!1,n)}set(n,a,e,t){let l=n[a];if(!this._isShallow){const c=ya(l);if(!Sn(e)&&!ya(e)&&(l=ws(l),e=ws(e)),!rs(n)&&zs(l)&&!zs(e))return c||(l.value=e),!0}const o=rs(n)&&po(a)?Number(a)<n.length:Ss(n,a),p=Reflect.set(n,a,e,zs(n)?n:t);return n===ws(t)&&(o?ba(e,l)&&Yn(n,"set",a,e):Yn(n,"add",a,e)),p}deleteProperty(n,a){const e=Ss(n,a);n[a];const t=Reflect.deleteProperty(n,a);return t&&e&&Yn(n,"delete",a,void 0),t}has(n,a){const e=Reflect.has(n,a);return(!la(a)||!hr.has(a))&&pn(n,"has",a),e}ownKeys(n){return pn(n,"iterate",rs(n)?"length":Ma),Reflect.ownKeys(n)}}class ah extends dr{constructor(n=!1){super(!0,n)}set(n,a){return!0}deleteProperty(n,a){return!0}}const eh=new mr,th=new ah,lh=new mr(!0);const kl=s=>s,Xe=s=>Reflect.getPrototypeOf(s);function oh(s,n,a){return function(...e){const t=this.__v_raw,l=ws(t),o=Ka(l),p=s==="entries"||s===Symbol.iterator&&o,c=s==="keys"&&o,r=t[s](...e),i=a?kl:n?jt:Zs;return!n&&pn(l,"iterate",c?Cl:Ma),{next(){const{value:u,done:h}=r.next();return h?{value:u,done:h}:{value:p?[i(u[0]),i(u[1])]:i(u),done:h}},[Symbol.iterator](){return this}}}}function Ye(s){return function(...n){return s==="delete"?!1:s==="clear"?void 0:this}}function ph(s,n){const a={get(t){const l=this.__v_raw,o=ws(l),p=ws(t);s||(ba(t,p)&&pn(o,"get",t),pn(o,"get",p));const{has:c}=Xe(o),r=n?kl:s?jt:Zs;if(c.call(o,t))return r(l.get(t));if(c.call(o,p))return r(l.get(p));l!==o&&l.get(t)},get size(){const t=this.__v_raw;return!s&&pn(ws(t),"iterate",Ma),t.size},has(t){const l=this.__v_raw,o=ws(l),p=ws(t);return s||(ba(t,p)&&pn(o,"has",t),pn(o,"has",p)),t===p?l.has(t):l.has(t)||l.has(p)},forEach(t,l){const o=this,p=o.__v_raw,c=ws(p),r=n?kl:s?jt:Zs;return!s&&pn(c,"iterate",Ma),p.forEach((i,u)=>t.call(l,r(i),r(u),o))}};return Ys(a,s?{add:Ye("add"),set:Ye("set"),delete:Ye("delete"),clear:Ye("clear")}:{add(t){!n&&!Sn(t)&&!ya(t)&&(t=ws(t));const l=ws(this);return Xe(l).has.call(l,t)||(l.add(t),Yn(l,"add",t,t)),this},set(t,l){!n&&!Sn(l)&&!ya(l)&&(l=ws(l));const o=ws(this),{has:p,get:c}=Xe(o);let r=p.call(o,t);r||(t=ws(t),r=p.call(o,t));const i=c.call(o,t);return o.set(t,l),r?ba(l,i)&&Yn(o,"set",t,l):Yn(o,"add",t,l),this},delete(t){const l=ws(this),{has:o,get:p}=Xe(l);let c=o.call(l,t);c||(t=ws(t),c=o.call(l,t)),p&&p.call(l,t);const r=l.delete(t);return c&&Yn(l,"delete",t,void 0),r},clear(){const t=ws(this),l=t.size!==0,o=t.clear();return l&&Yn(t,"clear",void 0,void 0),o}}),["keys","values","entries",Symbol.iterator].forEach(t=>{a[t]=oh(t,s,n)}),a}function fo(s,n){const a=ph(s,n);return(e,t,l)=>t==="__v_isReactive"?!s:t==="__v_isReadonly"?s:t==="__v_raw"?e:Reflect.get(Ss(a,t)&&t in e?a:e,t,l)}const ch={get:fo(!1,!1)},rh={get:fo(!1,!0)},ih={get:fo(!0,!1)};const jr=new WeakMap,fr=new WeakMap,gr=new WeakMap,uh=new WeakMap;function hh(s){switch(s){case"Object":case"Array":return 1;case"Map":case"Set":case"WeakMap":case"WeakSet":return 2;default:return 0}}function dh(s){return s.__v_skip||!Object.isExtensible(s)?0:hh($u(s))}function $e(s){return ya(s)?s:go(s,!1,eh,ch,jr)}function br(s){return go(s,!1,lh,rh,fr)}function xl(s){return go(s,!0,th,ih,gr)}function go(s,n,a,e,t){if(!Ms(s)||s.__v_raw&&!(n&&s.__v_isReactive))return s;const l=dh(s);if(l===0)return s;const o=t.get(s);if(o)return o;const p=new Proxy(s,l===2?e:a);return t.set(s,p),p}function _a(s){return ya(s)?_a(s.__v_raw):!!(s&&s.__v_isReactive)}function ya(s){return!!(s&&s.__v_isReadonly)}function Sn(s){return!!(s&&s.__v_isShallow)}function bo(s){return s?!!s.__v_raw:!1}function ws(s){const n=s&&s.__v_raw;return n?ws(n):s}function _o(s){return!Ss(s,"__v_skip")&&Object.isExtensible(s)&&Jc(s,"__v_skip",!0),s}const Zs=s=>Ms(s)?$e(s):s,jt=s=>Ms(s)?xl(s):s;function zs(s){return s?s.__v_isRef===!0:!1}function ps(s){return _r(s,!1)}function yo(s){return _r(s,!0)}function _r(s,n){return zs(s)?s:new mh(s,n)}class mh{constructor(n,a){this.dep=new jo,this.__v_isRef=!0,this.__v_isShallow=!1,this._rawValue=a?n:ws(n),this._value=a?n:Zs(n),this.__v_isShallow=a}get value(){return this.dep.track(),this._value}set value(n){const a=this._rawValue,e=this.__v_isShallow||Sn(n)||ya(n);n=e?n:ws(n),ba(n,a)&&(this._rawValue=n,this._value=e?n:Zs(n),this.dep.trigger())}}function Y(s){return zs(s)?s.value:s}function jh(s){return ds(s)?s():Y(s)}const fh={get:(s,n,a)=>n==="__v_raw"?s:Y(Reflect.get(s,n,a)),set:(s,n,a,e)=>{const t=s[n];return zs(t)&&!zs(a)?(t.value=a,!0):Reflect.set(s,n,a,e)}};function yr(s){return _a(s)?s:new Proxy(s,fh)}function gh(s){const n=rs(s)?new Array(s.length):{};for(const a in s)n[a]=_h(s,a);return n}class bh{constructor(n,a,e){this._object=n,this._key=a,this._defaultValue=e,this.__v_isRef=!0,this._value=void 0}get value(){const n=this._object[this._key];return this._value=n===void 0?this._defaultValue:n}set value(n){this._object[this._key]=n}get dep(){return Qu(ws(this._object),this._key)}}function _h(s,n,a){const e=s[n];return zs(e)?e:new bh(s,n,a)}class yh{constructor(n,a,e){this.fn=n,this.setter=a,this._value=void 0,this.dep=new jo(this),this.__v_isRef=!0,this.deps=void 0,this.depsTail=void 0,this.flags=16,this.globalVersion=ke-1,this.next=void 0,this.effect=this,this.__v_isReadonly=!a,this.isSSR=e}notify(){if(this.flags|=16,!(this.flags&8)&&Ls!==this)return or(this,!0),!0}get value(){const n=this.dep.track();return rr(this),n&&(n.version=this.dep.version),this._value}set value(n){this.setter&&this.setter(n)}}function vh(s,n,a=!1){let e,t;return ds(s)?e=s:(e=s.get,t=s.set),new yh(e,t,a)}const Qe={},ft=new WeakMap;let Ra;function wh(s,n=!1,a=Ra){if(a){let e=ft.get(a);e||ft.set(a,e=[]),e.push(s)}}function Ch(s,n,a=As){const{immediate:e,deep:t,once:l,scheduler:o,augmentJob:p,call:c}=a,r=w=>t?w:Sn(w)||t===!1||t===0?Qn(w,1):Qn(w);let i,u,h,d,f=!1,m=!1;if(zs(s)?(u=()=>s.value,f=Sn(s)):_a(s)?(u=()=>r(s),f=!0):rs(s)?(m=!0,f=s.some(w=>_a(w)||Sn(w)),u=()=>s.map(w=>{if(zs(w))return w.value;if(_a(w))return r(w);if(ds(w))return c?c(w,2):w()})):ds(s)?n?u=c?()=>c(s,2):s:u=()=>{if(h){sa();try{h()}finally{na()}}const w=Ra;Ra=i;try{return c?c(s,3,[d]):s(d)}finally{Ra=w}}:u=qn,n&&t){const w=u,R=t===!0?1/0:t;u=()=>Qn(w(),R)}const k=er(),j=()=>{i.stop(),k&&k.active&&oo(k.effects,i)};if(l&&n){const w=n;n=(...R)=>{w(...R),j()}}let v=m?new Array(s.length).fill(Qe):Qe;const y=w=>{if(!(!(i.flags&1)||!i.dirty&&!w))if(n){const R=i.run();if(t||f||(m?R.some((N,S)=>ba(N,v[S])):ba(R,v))){h&&h();const N=Ra;Ra=i;try{const S=[R,v===Qe?void 0:m&&v[0]===Qe?[]:v,d];v=R,c?c(n,3,S):n(...S)}finally{Ra=N}}}else i.run()};return p&&p(y),i=new tr(u),i.scheduler=o?()=>o(y,!1):y,d=w=>wh(w,!1,i),h=i.onStop=()=>{const w=ft.get(i);if(w){if(c)c(w,4);else for(const R of w)R();ft.delete(i)}},n?e?y(!0):v=i.run():o?o(y.bind(null,!0),!0):i.run(),j.pause=i.pause.bind(i),j.resume=i.resume.bind(i),j.stop=j,j}function Qn(s,n=1/0,a){if(n<=0||!Ms(s)||s.__v_skip||(a=a||new Map,(a.get(s)||0)>=n))return s;if(a.set(s,n),n--,zs(s))Qn(s.value,n,a);else if(rs(s))for(let e=0;e<s.length;e++)Qn(s[e],n,a);else if(Kc(s)||Ka(s))s.forEach(e=>{Qn(e,n,a)});else if(Qc(s)){for(const e in s)Qn(s[e],n,a);for(const e of Object.getOwnPropertySymbols(s))Object.prototype.propertyIsEnumerable.call(s,e)&&Qn(s[e],n,a)}return s}/**
* @vue/runtime-core v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/function Be(s,n,a,e){try{return e?s(...e):s()}catch(t){qe(t,n,a)}}function Ln(s,n,a,e){if(ds(s)){const t=Be(s,n,a,e);return t&&Xc(t)&&t.catch(l=>{qe(l,n,a)}),t}if(rs(s)){const t=[];for(let l=0;l<s.length;l++)t.push(Ln(s[l],n,a,e));return t}}function qe(s,n,a,e=!0){const t=n?n.vnode:null,{errorHandler:l,throwUnhandledErrorInProduction:o}=n&&n.appContext.config||As;if(n){let p=n.parent;const c=n.proxy,r=`https://vuejs.org/error-reference/#runtime-${a}`;for(;p;){const i=p.ec;if(i){for(let u=0;u<i.length;u++)if(i[u](s,c,r)===!1)return}p=p.parent}if(l){sa(),Be(l,null,10,[s,c,r]),na();return}}kh(s,a,t,e,o)}function kh(s,n,a,e=!0,t=!1){if(t)throw s;console.error(s)}const fn=[];let Fn=-1;const Xa=[];let ma=null,Ua=0;const vr=Promise.resolve();let gt=null;function vn(s){const n=gt||vr;return s?n.then(this?s.bind(this):s):n}function xh(s){let n=Fn+1,a=fn.length;for(;n<a;){const e=n+a>>>1,t=fn[e],l=Pe(t);l<s||l===s&&t.flags&2?n=e+1:a=e}return n}function vo(s){if(!(s.flags&1)){const n=Pe(s),a=fn[fn.length-1];!a||!(s.flags&2)&&n>=Pe(a)?fn.push(s):fn.splice(xh(n),0,s),s.flags|=1,wr()}}function wr(){gt||(gt=vr.then(kr))}function Ph(s){rs(s)?Xa.push(...s):ma&&s.id===-1?ma.splice(Ua+1,0,s):s.flags&1||(Xa.push(s),s.flags|=1),wr()}function Jo(s,n,a=Fn+1){for(;a<fn.length;a++){const e=fn[a];if(e&&e.flags&2){if(s&&e.id!==s.uid)continue;fn.splice(a,1),a--,e.flags&4&&(e.flags&=-2),e(),e.flags&4||(e.flags&=-2)}}}function Cr(s){if(Xa.length){const n=[...new Set(Xa)].sort((a,e)=>Pe(a)-Pe(e));if(Xa.length=0,ma){ma.push(...n);return}for(ma=n,Ua=0;Ua<ma.length;Ua++){const a=ma[Ua];a.flags&4&&(a.flags&=-2),a.flags&8||a(),a.flags&=-2}ma=null,Ua=0}}const Pe=s=>s.id==null?s.flags&2?-1:1/0:s.id;function kr(s){try{for(Fn=0;Fn<fn.length;Fn++){const n=fn[Fn];n&&!(n.flags&8)&&(n.flags&4&&(n.flags&=-2),Be(n,n.i,n.i?15:14),n.flags&4||(n.flags&=-2))}}finally{for(;Fn<fn.length;Fn++){const n=fn[Fn];n&&(n.flags&=-2)}Fn=-1,fn.length=0,Cr(),gt=null,(fn.length||Xa.length)&&kr()}}let nn=null,xr=null;function bt(s){const n=nn;return nn=s,xr=s&&s.type.__scopeId||null,n}function an(s,n=nn,a){if(!n||s._n)return s;const e=(...t)=>{e._d&&vt(-1);const l=bt(n);let o;try{o=s(...t)}finally{bt(l),e._d&&vt(1)}return o};return e._n=!0,e._c=!0,e._d=!0,e}function gn(s,n){if(nn===null)return s;const a=$t(nn),e=s.dirs||(s.dirs=[]);for(let t=0;t<n.length;t++){let[l,o,p,c=As]=n[t];l&&(ds(l)&&(l={mounted:l,updated:l}),l.deep&&Qn(o),e.push({dir:l,instance:a,value:o,oldValue:void 0,arg:p,modifiers:c}))}return s}function Sa(s,n,a,e){const t=s.dirs,l=n&&n.dirs;for(let o=0;o<t.length;o++){const p=t[o];l&&(p.oldValue=l[o].value);let c=p.dir[e];c&&(sa(),Ln(c,a,8,[s.el,p,s,n]),na())}}const Pr=Symbol("_vte"),Sr=s=>s.__isTeleport,fe=s=>s&&(s.disabled||s.disabled===""),Zo=s=>s&&(s.defer||s.defer===""),sp=s=>typeof SVGElement<"u"&&s instanceof SVGElement,np=s=>typeof MathMLElement=="function"&&s instanceof MathMLElement,Pl=(s,n)=>{const a=s&&s.to;return Us(a)?n?n(a):null:a},Er={name:"Teleport",__isTeleport:!0,process(s,n,a,e,t,l,o,p,c,r){const{mc:i,pc:u,pbc:h,o:{insert:d,querySelector:f,createText:m,createComment:k}}=r,j=fe(n.props);let{shapeFlag:v,children:y,dynamicChildren:w}=n;if(s==null){const R=n.el=m(""),N=n.anchor=m("");d(R,a,e),d(N,a,e);const S=(T,z)=>{v&16&&i(y,T,z,t,l,o,p,c)},A=()=>{const T=n.target=Pl(n.props,f),z=Tr(T,n,m,d);T&&(o!=="svg"&&sp(T)?o="svg":o!=="mathml"&&np(T)&&(o="mathml"),t&&t.isCE&&(t.ce._teleportTargets||(t.ce._teleportTargets=new Set)).add(T),j||(S(T,z),pt(n,!1)))};j&&(S(a,N),pt(n,!0)),Zo(n.props)?(n.el.__isMounted=!1,mn(()=>{A(),delete n.el.__isMounted},l)):A()}else{if(Zo(n.props)&&s.el.__isMounted===!1){mn(()=>{Er.process(s,n,a,e,t,l,o,p,c,r)},l);return}n.el=s.el,n.targetStart=s.targetStart;const R=n.anchor=s.anchor,N=n.target=s.target,S=n.targetAnchor=s.targetAnchor,A=fe(s.props),T=A?a:N,z=A?R:S;if(o==="svg"||sp(N)?o="svg":(o==="mathml"||np(N))&&(o="mathml"),w?(h(s.dynamicChildren,w,T,t,l,o,p),Eo(s,n,!0)):c||u(s,n,T,z,t,l,o,p,!1),j)A?n.props&&s.props&&n.props.to!==s.props.to&&(n.props.to=s.props.to):Je(n,a,R,r,1);else if((n.props&&n.props.to)!==(s.props&&s.props.to)){const I=n.target=Pl(n.props,f);I&&Je(n,I,null,r,0)}else A&&Je(n,N,S,r,1);pt(n,j)}},remove(s,n,a,{um:e,o:{remove:t}},l){const{shapeFlag:o,children:p,anchor:c,targetStart:r,targetAnchor:i,target:u,props:h}=s;if(u&&(t(r),t(i)),l&&t(c),o&16){const d=l||!fe(h);for(let f=0;f<p.length;f++){const m=p[f];e(m,n,a,d,!!m.dynamicChildren)}}},move:Je,hydrate:Sh};function Je(s,n,a,{o:{insert:e},m:t},l=2){l===0&&e(s.targetAnchor,n,a);const{el:o,anchor:p,shapeFlag:c,children:r,props:i}=s,u=l===2;if(u&&e(o,n,a),(!u||fe(i))&&c&16)for(let h=0;h<r.length;h++)t(r[h],n,a,2);u&&e(p,n,a)}function Sh(s,n,a,e,t,l,{o:{nextSibling:o,parentNode:p,querySelector:c,insert:r,createText:i}},u){function h(m,k,j,v){k.anchor=u(o(m),k,p(m),a,e,t,l),k.targetStart=j,k.targetAnchor=v}const d=n.target=Pl(n.props,c),f=fe(n.props);if(d){const m=d._lpa||d.firstChild;if(n.shapeFlag&16)if(f)h(s,n,m,m&&o(m));else{n.anchor=o(s);let k=m;for(;k;){if(k&&k.nodeType===8){if(k.data==="teleport start anchor")n.targetStart=k;else if(k.data==="teleport anchor"){n.targetAnchor=k,d._lpa=n.targetAnchor&&o(n.targetAnchor);break}}k=o(k)}n.targetAnchor||Tr(d,n,i,r),u(m&&o(m),n,d,a,e,t,l)}pt(n,f)}else f&&n.shapeFlag&16&&h(s,n,s,o(s));return n.anchor&&o(n.anchor)}const C2=Er;function pt(s,n){const a=s.ctx;if(a&&a.ut){let e,t;for(n?(e=s.el,t=s.anchor):(e=s.targetStart,t=s.targetAnchor);e&&e!==t;)e.nodeType===1&&e.setAttribute("data-v-owner",a.uid),e=e.nextSibling;a.ut()}}function Tr(s,n,a,e){const t=n.targetStart=a(""),l=n.targetAnchor=a("");return t[Pr]=l,s&&(e(t,s),e(l,s)),l}const Xn=Symbol("_leaveCb"),Ze=Symbol("_enterCb");function Eh(){const s={isMounted:!1,isLeaving:!1,isUnmounting:!1,leavingVNodes:new Map};return kn(()=>{s.isMounted=!0}),Dn(()=>{s.isUnmounting=!0}),s}const xn=[Function,Array],Ar={mode:String,appear:Boolean,persisted:Boolean,onBeforeEnter:xn,onEnter:xn,onAfterEnter:xn,onEnterCancelled:xn,onBeforeLeave:xn,onLeave:xn,onAfterLeave:xn,onLeaveCancelled:xn,onBeforeAppear:xn,onAppear:xn,onAfterAppear:xn,onAppearCancelled:xn},Rr=s=>{const n=s.subTree;return n.component?Rr(n.component):n},Th={name:"BaseTransition",props:Ar,setup(s,{slots:n}){const a=ca(),e=Eh();return()=>{const t=n.default&&Dr(n.default(),!0);if(!t||!t.length)return;const l=Lr(t),o=ws(s),{mode:p}=o;if(e.isLeaving)return al(l);const c=ap(l);if(!c)return al(l);let r=Sl(c,o,e,a,u=>r=u);c.type!==cn&&Se(c,r);let i=a.subTree&&ap(a.subTree);if(i&&i.type!==cn&&!La(i,c)&&Rr(a).type!==cn){let u=Sl(i,o,e,a);if(Se(i,u),p==="out-in"&&c.type!==cn)return e.isLeaving=!0,u.afterLeave=()=>{e.isLeaving=!1,a.job.flags&8||a.update(),delete u.afterLeave,i=void 0},al(l);p==="in-out"&&c.type!==cn?u.delayLeave=(h,d,f)=>{const m=Mr(e,i);m[String(i.key)]=i,h[Xn]=()=>{d(),h[Xn]=void 0,delete r.delayedLeave,i=void 0},r.delayedLeave=()=>{f(),delete r.delayedLeave,i=void 0}}:i=void 0}else i&&(i=void 0);return l}}};function Lr(s){let n=s[0];if(s.length>1){for(const a of s)if(a.type!==cn){n=a;break}}return n}const Ah=Th;function Mr(s,n){const{leavingVNodes:a}=s;let e=a.get(n.type);return e||(e=Object.create(null),a.set(n.type,e)),e}function Sl(s,n,a,e,t){const{appear:l,mode:o,persisted:p=!1,onBeforeEnter:c,onEnter:r,onAfterEnter:i,onEnterCancelled:u,onBeforeLeave:h,onLeave:d,onAfterLeave:f,onLeaveCancelled:m,onBeforeAppear:k,onAppear:j,onAfterAppear:v,onAppearCancelled:y}=n,w=String(s.key),R=Mr(a,s),N=(T,z)=>{T&&Ln(T,e,9,z)},S=(T,z)=>{const I=z[1];N(T,z),rs(T)?T.every(D=>D.length<=1)&&I():T.length<=1&&I()},A={mode:o,persisted:p,beforeEnter(T){let z=c;if(!a.isMounted)if(l)z=k||c;else return;T[Xn]&&T[Xn](!0);const I=R[w];I&&La(s,I)&&I.el[Xn]&&I.el[Xn](),N(z,[T])},enter(T){let z=r,I=i,D=u;if(!a.isMounted)if(l)z=j||r,I=v||i,D=y||u;else return;let W=!1;const ls=T[Ze]=ms=>{W||(W=!0,ms?N(D,[T]):N(I,[T]),A.delayedLeave&&A.delayedLeave(),T[Ze]=void 0)};z?S(z,[T,ls]):ls()},leave(T,z){const I=String(s.key);if(T[Ze]&&T[Ze](!0),a.isUnmounting)return z();N(h,[T]);let D=!1;const W=T[Xn]=ls=>{D||(D=!0,z(),ls?N(m,[T]):N(f,[T]),T[Xn]=void 0,R[I]===s&&delete R[I])};R[I]=s,d?S(d,[T,W]):W()},clone(T){const z=Sl(T,n,a,e,t);return t&&t(z),z}};return A}function al(s){if(ze(s))return s=va(s),s.children=null,s}function ap(s){if(!ze(s))return Sr(s.type)&&s.children?Lr(s.children):s;if(s.component)return s.component.subTree;const{shapeFlag:n,children:a}=s;if(a){if(n&16)return a[0];if(n&32&&ds(a.default))return a.default()}}function Se(s,n){s.shapeFlag&6&&s.component?(s.transition=n,Se(s.component.subTree,n)):s.shapeFlag&128?(s.ssContent.transition=n.clone(s.ssContent),s.ssFallback.transition=n.clone(s.ssFallback)):s.transition=n}function Dr(s,n=!1,a){let e=[],t=0;for(let l=0;l<s.length;l++){let o=s[l];const p=a==null?o.key:String(a)+String(o.key!=null?o.key:l);o.type===js?(o.patchFlag&128&&t++,e=e.concat(Dr(o.children,n,p))):(n||o.type!==cn)&&e.push(p!=null?va(o,{key:p}):o)}if(t>1)for(let l=0;l<e.length;l++)e[l].patchFlag=-2;return e}function Es(s,n){return ds(s)?Ys({name:s.name},n,{setup:s}):s}function wo(s){s.ids=[s.ids[0]+s.ids[2]+++"-",0,0]}function ct(s){const n=ca(),a=yo(null);if(n){const t=n.refs===As?n.refs={}:n.refs;Object.defineProperty(t,s,{enumerable:!0,get:()=>a.value,set:l=>a.value=l})}return a}const _t=new WeakMap;function ge(s,n,a,e,t=!1){if(rs(s)){s.forEach((f,m)=>ge(f,n&&(rs(n)?n[m]:n),a,e,t));return}if(Ya(e)&&!t){e.shapeFlag&512&&e.type.__asyncResolved&&e.component.subTree.component&&ge(s,n,a,e.component.subTree);return}const l=e.shapeFlag&4?$t(e.component):e.el,o=t?null:l,{i:p,r:c}=s,r=n&&n.r,i=p.refs===As?p.refs={}:p.refs,u=p.setupState,h=ws(u),d=u===As?Wc:f=>Ss(h,f);if(r!=null&&r!==c){if(ep(n),Us(r))i[r]=null,d(r)&&(u[r]=null);else if(zs(r)){r.value=null;const f=n;f.k&&(i[f.k]=null)}}if(ds(c))Be(c,p,12,[o,i]);else{const f=Us(c),m=zs(c);if(f||m){const k=()=>{if(s.f){const j=f?d(c)?u[c]:i[c]:c.value;if(t)rs(j)&&oo(j,l);else if(rs(j))j.includes(l)||j.push(l);else if(f)i[c]=[l],d(c)&&(u[c]=i[c]);else{const v=[l];c.value=v,s.k&&(i[s.k]=v)}}else f?(i[c]=o,d(c)&&(u[c]=o)):m&&(c.value=o,s.k&&(i[s.k]=o))};if(o){const j=()=>{k(),_t.delete(s)};j.id=-1,_t.set(s,j),mn(j,a)}else ep(s),k()}}}function ep(s){const n=_t.get(s);n&&(n.flags|=8,_t.delete(s))}const tp=s=>s.nodeType===8;Ot().requestIdleCallback;Ot().cancelIdleCallback;function Rh(s,n){if(tp(s)&&s.data==="["){let a=1,e=s.nextSibling;for(;e;){if(e.nodeType===1){if(n(e)===!1)break}else if(tp(e))if(e.data==="]"){if(--a===0)break}else e.data==="["&&a++;e=e.nextSibling}}else n(s)}const Ya=s=>!!s.type.__asyncLoader;function Lh(s){ds(s)&&(s={loader:s});const{loader:n,loadingComponent:a,errorComponent:e,delay:t=200,hydrate:l,timeout:o,suspensible:p=!0,onError:c}=s;let r=null,i,u=0;const h=()=>(u++,r=null,d()),d=()=>{let f;return r||(f=r=n().catch(m=>{if(m=m instanceof Error?m:new Error(String(m)),c)return new Promise((k,j)=>{c(m,()=>k(h()),()=>j(m),u+1)});throw m}).then(m=>f!==r&&r?r:(m&&(m.__esModule||m[Symbol.toStringTag]==="Module")&&(m=m.default),i=m,m)))};return Es({name:"AsyncComponentWrapper",__asyncLoader:d,__asyncHydrate(f,m,k){let j=!1;(m.bu||(m.bu=[])).push(()=>j=!0);const v=()=>{j||k()},y=l?()=>{const w=l(v,R=>Rh(f,R));w&&(m.bum||(m.bum=[])).push(w)}:v;i?y():d().then(()=>!m.isUnmounted&&y())},get __asyncResolved(){return i},setup(){const f=sn;if(wo(f),i)return()=>st(i,f);const m=y=>{r=null,qe(y,f,13,!e)};if(p&&f.suspense||Ja)return d().then(y=>()=>st(y,f)).catch(y=>(m(y),()=>e?fs(e,{error:y}):null));const k=ps(!1),j=ps(),v=ps(!!t);return t&&setTimeout(()=>{v.value=!1},t),o!=null&&setTimeout(()=>{if(!k.value&&!j.value){const y=new Error(`Async component timed out after ${o}ms.`);m(y),j.value=y}},o),d().then(()=>{k.value=!0,f.parent&&ze(f.parent.vnode)&&f.parent.update()}).catch(y=>{m(y),j.value=y}),()=>{if(k.value&&i)return st(i,f);if(j.value&&e)return fs(e,{error:j.value});if(a&&!v.value)return st(a,f)}}})}function st(s,n){const{ref:a,props:e,children:t,ce:l}=n.vnode,o=fs(s,e,t);return o.ref=a,o.ce=l,delete n.vnode.ce,o}const ze=s=>s.type.__isKeepAlive;function Or(s,n){Nr(s,"a",n)}function Ir(s,n){Nr(s,"da",n)}function Nr(s,n,a=sn){const e=s.__wdc||(s.__wdc=()=>{let t=a;for(;t;){if(t.isDeactivated)return;t=t.parent}return s()});if(Nt(n,e,a),a){let t=a.parent;for(;t&&t.parent;)ze(t.parent.vnode)&&Mh(e,n,a,t),t=t.parent}}function Mh(s,n,a,e){const t=Nt(n,s,e,!0);Co(()=>{oo(e[n],t)},a)}function Nt(s,n,a=sn,e=!1){if(a){const t=a[s]||(a[s]=[]),l=n.__weh||(n.__weh=(...o)=>{sa();const p=He(a),c=Ln(n,a,s,o);return p(),na(),c});return e?t.unshift(l):t.push(l),l}}const pa=s=>(n,a=sn)=>{(!Ja||s==="sp")&&Nt(s,(...e)=>n(...e),a)},Dh=pa("bm"),kn=pa("m"),Oh=pa("bu"),Ih=pa("u"),Dn=pa("bum"),Co=pa("um"),Nh=pa("sp"),Fh=pa("rtg"),$h=pa("rtc");function Bh(s,n=sn){Nt("ec",s,n)}const ko="components",qh="directives";function aa(s,n){return xo(ko,s,!0,n)||s}const Fr=Symbol.for("v-ndc");function zh(s){return Us(s)?xo(ko,s,!1)||s:s||Fr}function ae(s){return xo(qh,s)}function xo(s,n,a=!0,e=!1){const t=nn||sn;if(t){const l=t.type;if(s===ko){const p=Ad(l,!1);if(p&&(p===n||p===Tn(n)||p===Dt(Tn(n))))return l}const o=lp(t[s]||l[s],n)||lp(t.appContext[s],n);return!o&&e?l:o}}function lp(s,n){return s&&(s[n]||s[Tn(n)]||s[Dt(Tn(n))])}function Is(s,n,a,e){let t;const l=a,o=rs(s);if(o||Us(s)){const p=o&&_a(s);let c=!1,r=!1;p&&(c=!Sn(s),r=ya(s),s=It(s)),t=new Array(s.length);for(let i=0,u=s.length;i<u;i++)t[i]=n(c?r?jt(Zs(s[i])):Zs(s[i]):s[i],i,void 0,l)}else if(typeof s=="number"){t=new Array(s);for(let p=0;p<s;p++)t[p]=n(p+1,p,void 0,l)}else if(Ms(s))if(s[Symbol.iterator])t=Array.from(s,(p,c)=>n(p,c,void 0,l));else{const p=Object.keys(s);t=new Array(p.length);for(let c=0,r=p.length;c<r;c++){const i=p[c];t[c]=n(s[i],i,c,l)}}else t=[];return t}function $r(s,n,a={},e,t){if(nn.ce||nn.parent&&Ya(nn.parent)&&nn.parent.ce){const r=Object.keys(a).length>0;return O(),rn(js,null,[fs("slot",a,e)],r?-2:64)}let l=s[n];l&&l._c&&(l._d=!1),O();const o=l&&Br(l(a)),p=a.key||o&&o.key,c=rn(js,{key:(p&&!la(p)?p:`_${n}`)+(!o&&e?"_fb":"")},o||[],o&&s._===1?64:-2);return l&&l._c&&(l._d=!0),c}function Br(s){return s.some(n=>Te(n)?!(n.type===cn||n.type===js&&!Br(n.children)):!0)?s:null}const El=s=>s?li(s)?$t(s):El(s.parent):null,be=Ys(Object.create(null),{$:s=>s,$el:s=>s.vnode.el,$data:s=>s.data,$props:s=>s.props,$attrs:s=>s.attrs,$slots:s=>s.slots,$refs:s=>s.refs,$parent:s=>El(s.parent),$root:s=>El(s.root),$host:s=>s.ce,$emit:s=>s.emit,$options:s=>zr(s),$forceUpdate:s=>s.f||(s.f=()=>{vo(s.update)}),$nextTick:s=>s.n||(s.n=vn.bind(s.proxy)),$watch:s=>id.bind(s)}),el=(s,n)=>s!==As&&!s.__isScriptSetup&&Ss(s,n),Uh={get({_:s},n){if(n==="__v_skip")return!0;const{ctx:a,setupState:e,data:t,props:l,accessCache:o,type:p,appContext:c}=s;let r;if(n[0]!=="$"){const d=o[n];if(d!==void 0)switch(d){case 1:return e[n];case 2:return t[n];case 4:return a[n];case 3:return l[n]}else{if(el(e,n))return o[n]=1,e[n];if(t!==As&&Ss(t,n))return o[n]=2,t[n];if((r=s.propsOptions[0])&&Ss(r,n))return o[n]=3,l[n];if(a!==As&&Ss(a,n))return o[n]=4,a[n];Tl&&(o[n]=0)}}const i=be[n];let u,h;if(i)return n==="$attrs"&&pn(s.attrs,"get",""),i(s);if((u=p.__cssModules)&&(u=u[n]))return u;if(a!==As&&Ss(a,n))return o[n]=4,a[n];if(h=c.config.globalProperties,Ss(h,n))return h[n]},set({_:s},n,a){const{data:e,setupState:t,ctx:l}=s;return el(t,n)?(t[n]=a,!0):e!==As&&Ss(e,n)?(e[n]=a,!0):Ss(s.props,n)||n[0]==="$"&&n.slice(1)in s?!1:(l[n]=a,!0)},has({_:{data:s,setupState:n,accessCache:a,ctx:e,appContext:t,propsOptions:l,type:o}},p){let c,r;return!!(a[p]||s!==As&&p[0]!=="$"&&Ss(s,p)||el(n,p)||(c=l[0])&&Ss(c,p)||Ss(e,p)||Ss(be,p)||Ss(t.config.globalProperties,p)||(r=o.__cssModules)&&r[p])},defineProperty(s,n,a){return a.get!=null?s._.accessCache[n]=0:Ss(a,"value")&&this.set(s,n,a.value,null),Reflect.defineProperty(s,n,a)}};function op(s){return rs(s)?s.reduce((n,a)=>(n[a]=null,n),{}):s}let Tl=!0;function Hh(s){const n=zr(s),a=s.proxy,e=s.ctx;Tl=!1,n.beforeCreate&&pp(n.beforeCreate,s,"bc");const{data:t,computed:l,methods:o,watch:p,provide:c,inject:r,created:i,beforeMount:u,mounted:h,beforeUpdate:d,updated:f,activated:m,deactivated:k,beforeDestroy:j,beforeUnmount:v,destroyed:y,unmounted:w,render:R,renderTracked:N,renderTriggered:S,errorCaptured:A,serverPrefetch:T,expose:z,inheritAttrs:I,components:D,directives:W,filters:ls}=n;if(r&&Vh(r,e,null),o)for(const es in o){const cs=o[es];ds(cs)&&(e[es]=cs.bind(a))}if(t){const es=t.call(a,a);Ms(es)&&(s.data=$e(es))}if(Tl=!0,l)for(const es in l){const cs=l[es],$s=ds(cs)?cs.bind(a,a):ds(cs.get)?cs.get.bind(a,a):qn,Bs=!ds(cs)&&ds(cs.set)?cs.set.bind(a):qn,qs=J({get:$s,set:Bs});Object.defineProperty(e,es,{enumerable:!0,configurable:!0,get:()=>qs.value,set:os=>qs.value=os})}if(p)for(const es in p)qr(p[es],e,a,es);if(c){const es=ds(c)?c.call(a):c;Reflect.ownKeys(es).forEach(cs=>{rt(cs,es[cs])})}i&&pp(i,s,"c");function us(es,cs){rs(cs)?cs.forEach($s=>es($s.bind(a))):cs&&es(cs.bind(a))}if(us(Dh,u),us(kn,h),us(Oh,d),us(Ih,f),us(Or,m),us(Ir,k),us(Bh,A),us($h,N),us(Fh,S),us(Dn,v),us(Co,w),us(Nh,T),rs(z))if(z.length){const es=s.exposed||(s.exposed={});z.forEach(cs=>{Object.defineProperty(es,cs,{get:()=>a[cs],set:$s=>a[cs]=$s,enumerable:!0})})}else s.exposed||(s.exposed={});R&&s.render===qn&&(s.render=R),I!=null&&(s.inheritAttrs=I),D&&(s.components=D),W&&(s.directives=W),T&&wo(s)}function Vh(s,n,a=qn){rs(s)&&(s=Al(s));for(const e in s){const t=s[e];let l;Ms(t)?"default"in t?l=_n(t.from||e,t.default,!0):l=_n(t.from||e):l=_n(t),zs(l)?Object.defineProperty(n,e,{enumerable:!0,configurable:!0,get:()=>l.value,set:o=>l.value=o}):n[e]=l}}function pp(s,n,a){Ln(rs(s)?s.map(e=>e.bind(n.proxy)):s.bind(n.proxy),n,a)}function qr(s,n,a,e){let t=e.includes(".")?si(a,e):()=>a[e];if(Us(s)){const l=n[s];ds(l)&&Ws(t,l)}else if(ds(s))Ws(t,s.bind(a));else if(Ms(s))if(rs(s))s.forEach(l=>qr(l,n,a,e));else{const l=ds(s.handler)?s.handler.bind(a):n[s.handler];ds(l)&&Ws(t,l,s)}}function zr(s){const n=s.type,{mixins:a,extends:e}=n,{mixins:t,optionsCache:l,config:{optionMergeStrategies:o}}=s.appContext,p=l.get(n);let c;return p?c=p:!t.length&&!a&&!e?c=n:(c={},t.length&&t.forEach(r=>yt(c,r,o,!0)),yt(c,n,o)),Ms(n)&&l.set(n,c),c}function yt(s,n,a,e=!1){const{mixins:t,extends:l}=n;l&&yt(s,l,a,!0),t&&t.forEach(o=>yt(s,o,a,!0));for(const o in n)if(!(e&&o==="expose")){const p=Gh[o]||a&&a[o];s[o]=p?p(s[o],n[o]):n[o]}return s}const Gh={data:cp,props:rp,emits:rp,methods:he,computed:he,beforeCreate:hn,created:hn,beforeMount:hn,mounted:hn,beforeUpdate:hn,updated:hn,beforeDestroy:hn,beforeUnmount:hn,destroyed:hn,unmounted:hn,activated:hn,deactivated:hn,errorCaptured:hn,serverPrefetch:hn,components:he,directives:he,watch:Kh,provide:cp,inject:Wh};function cp(s,n){return n?s?function(){return Ys(ds(s)?s.call(this,this):s,ds(n)?n.call(this,this):n)}:n:s}function Wh(s,n){return he(Al(s),Al(n))}function Al(s){if(rs(s)){const n={};for(let a=0;a<s.length;a++)n[s[a]]=s[a];return n}return s}function hn(s,n){return s?[...new Set([].concat(s,n))]:n}function he(s,n){return s?Ys(Object.create(null),s,n):n}function rp(s,n){return s?rs(s)&&rs(n)?[...new Set([...s,...n])]:Ys(Object.create(null),op(s),op(n??{})):n}function Kh(s,n){if(!s)return n;if(!n)return s;const a=Ys(Object.create(null),s);for(const e in n)a[e]=hn(s[e],n[e]);return a}function Ur(){return{app:null,config:{isNativeTag:Wc,performance:!1,globalProperties:{},optionMergeStrategies:{},errorHandler:void 0,warnHandler:void 0,compilerOptions:{}},mixins:[],components:{},directives:{},provides:Object.create(null),optionsCache:new WeakMap,propsCache:new WeakMap,emitsCache:new WeakMap}}let Xh=0;function Yh(s,n){return function(e,t=null){ds(e)||(e=Ys({},e)),t!=null&&!Ms(t)&&(t=null);const l=Ur(),o=new WeakSet,p=[];let c=!1;const r=l.app={_uid:Xh++,_component:e,_props:t,_container:null,_context:l,_instance:null,version:Ld,get config(){return l.config},set config(i){},use(i,...u){return o.has(i)||(i&&ds(i.install)?(o.add(i),i.install(r,...u)):ds(i)&&(o.add(i),i(r,...u))),r},mixin(i){return l.mixins.includes(i)||l.mixins.push(i),r},component(i,u){return u?(l.components[i]=u,r):l.components[i]},directive(i,u){return u?(l.directives[i]=u,r):l.directives[i]},mount(i,u,h){if(!c){const d=r._ceVNode||fs(e,t);return d.appContext=l,h===!0?h="svg":h===!1&&(h=void 0),s(d,i,h),c=!0,r._container=i,i.__vue_app__=r,$t(d.component)}},onUnmount(i){p.push(i)},unmount(){c&&(Ln(p,r._instance,16),s(null,r._container),delete r._container.__vue_app__)},provide(i,u){return l.provides[i]=u,r},runWithContext(i){const u=Da;Da=r;try{return i()}finally{Da=u}}};return r}}let Da=null;function rt(s,n){if(sn){let a=sn.provides;const e=sn.parent&&sn.parent.provides;e===a&&(a=sn.provides=Object.create(e)),a[s]=n}}function _n(s,n,a=!1){const e=ca();if(e||Da){let t=Da?Da._context.provides:e?e.parent==null||e.ce?e.vnode.appContext&&e.vnode.appContext.provides:e.parent.provides:void 0;if(t&&s in t)return t[s];if(arguments.length>1)return a&&ds(n)?n.call(e&&e.proxy):n}}function Hr(){return!!(ca()||Da)}const Vr={},Gr=()=>Object.create(Vr),Wr=s=>Object.getPrototypeOf(s)===Vr;function Qh(s,n,a,e=!1){const t={},l=Gr();s.propsDefaults=Object.create(null),Kr(s,n,t,l);for(const o in s.propsOptions[0])o in t||(t[o]=void 0);a?s.props=e?t:br(t):s.type.props?s.props=t:s.props=l,s.attrs=l}function Jh(s,n,a,e){const{props:t,attrs:l,vnode:{patchFlag:o}}=s,p=ws(t),[c]=s.propsOptions;let r=!1;if((e||o>0)&&!(o&16)){if(o&8){const i=s.vnode.dynamicProps;for(let u=0;u<i.length;u++){let h=i[u];if(Ft(s.emitsOptions,h))continue;const d=n[h];if(c)if(Ss(l,h))d!==l[h]&&(l[h]=d,r=!0);else{const f=Tn(h);t[f]=Rl(c,p,f,d,s,!1)}else d!==l[h]&&(l[h]=d,r=!0)}}}else{Kr(s,n,t,l)&&(r=!0);let i;for(const u in p)(!n||!Ss(n,u)&&((i=Ca(u))===u||!Ss(n,i)))&&(c?a&&(a[u]!==void 0||a[i]!==void 0)&&(t[u]=Rl(c,p,u,void 0,s,!0)):delete t[u]);if(l!==p)for(const u in l)(!n||!Ss(n,u))&&(delete l[u],r=!0)}r&&Yn(s.attrs,"set","")}function Kr(s,n,a,e){const[t,l]=s.propsOptions;let o=!1,p;if(n)for(let c in n){if(de(c))continue;const r=n[c];let i;t&&Ss(t,i=Tn(c))?!l||!l.includes(i)?a[i]=r:(p||(p={}))[i]=r:Ft(s.emitsOptions,c)||(!(c in e)||r!==e[c])&&(e[c]=r,o=!0)}if(l){const c=ws(a),r=p||As;for(let i=0;i<l.length;i++){const u=l[i];a[u]=Rl(t,c,u,r[u],s,!Ss(r,u))}}return o}function Rl(s,n,a,e,t,l){const o=s[a];if(o!=null){const p=Ss(o,"default");if(p&&e===void 0){const c=o.default;if(o.type!==Function&&!o.skipFactory&&ds(c)){const{propsDefaults:r}=t;if(a in r)e=r[a];else{const i=He(t);e=r[a]=c.call(null,n),i()}}else e=c;t.ce&&t.ce._setProp(a,e)}o[0]&&(l&&!p?e=!1:o[1]&&(e===""||e===Ca(a))&&(e=!0))}return e}const Zh=new WeakMap;function Xr(s,n,a=!1){const e=a?Zh:n.propsCache,t=e.get(s);if(t)return t;const l=s.props,o={},p=[];let c=!1;if(!ds(s)){const i=u=>{c=!0;const[h,d]=Xr(u,n,!0);Ys(o,h),d&&p.push(...d)};!a&&n.mixins.length&&n.mixins.forEach(i),s.extends&&i(s.extends),s.mixins&&s.mixins.forEach(i)}if(!l&&!c)return Ms(s)&&e.set(s,Wa),Wa;if(rs(l))for(let i=0;i<l.length;i++){const u=Tn(l[i]);ip(u)&&(o[u]=As)}else if(l)for(const i in l){const u=Tn(i);if(ip(u)){const h=l[i],d=o[u]=rs(h)||ds(h)?{type:h}:Ys({},h),f=d.type;let m=!1,k=!0;if(rs(f))for(let j=0;j<f.length;++j){const v=f[j],y=ds(v)&&v.name;if(y==="Boolean"){m=!0;break}else y==="String"&&(k=!1)}else m=ds(f)&&f.name==="Boolean";d[0]=m,d[1]=k,(m||Ss(d,"default"))&&p.push(u)}}const r=[o,p];return Ms(s)&&e.set(s,r),r}function ip(s){return s[0]!=="$"&&!de(s)}const Po=s=>s==="_"||s==="_ctx"||s==="$stable",So=s=>rs(s)?s.map($n):[$n(s)],sd=(s,n,a)=>{if(n._n)return n;const e=an((...t)=>So(n(...t)),a);return e._c=!1,e},Yr=(s,n,a)=>{const e=s._ctx;for(const t in s){if(Po(t))continue;const l=s[t];if(ds(l))n[t]=sd(t,l,e);else if(l!=null){const o=So(l);n[t]=()=>o}}},Qr=(s,n)=>{const a=So(n);s.slots.default=()=>a},Jr=(s,n,a)=>{for(const e in n)(a||!Po(e))&&(s[e]=n[e])},nd=(s,n,a)=>{const e=s.slots=Gr();if(s.vnode.shapeFlag&32){const t=n._;t?(Jr(e,n,a),a&&Jc(e,"_",t,!0)):Yr(n,e)}else n&&Qr(s,n)},ad=(s,n,a)=>{const{vnode:e,slots:t}=s;let l=!0,o=As;if(e.shapeFlag&32){const p=n._;p?a&&p===1?l=!1:Jr(t,n,a):(l=!n.$stable,Yr(n,t)),o=n}else n&&(Qr(s,n),o={default:1});if(l)for(const p in t)!Po(p)&&o[p]==null&&delete t[p]},mn=bd;function ed(s){return td(s)}function td(s,n){const a=Ot();a.__VUE__=!0;const{insert:e,remove:t,patchProp:l,createElement:o,createText:p,createComment:c,setText:r,setElementText:i,parentNode:u,nextSibling:h,setScopeId:d=qn,insertStaticContent:f}=s,m=(_,C,E,B=null,H=null,U=null,g=void 0,b=null,P=!!C.dynamicChildren)=>{if(_===C)return;_&&!La(_,C)&&(B=F(_),os(_,H,U,!0),_=null),C.patchFlag===-2&&(P=!1,C.dynamicChildren=null);const{type:M,ref:Q,shapeFlag:G}=C;switch(M){case Ue:k(_,C,E,B);break;case cn:j(_,C,E,B);break;case ll:_==null&&v(C,E,B,g);break;case js:D(_,C,E,B,H,U,g,b,P);break;default:G&1?R(_,C,E,B,H,U,g,b,P):G&6?W(_,C,E,B,H,U,g,b,P):(G&64||G&128)&&M.process(_,C,E,B,H,U,g,b,P,as)}Q!=null&&H?ge(Q,_&&_.ref,U,C||_,!C):Q==null&&_&&_.ref!=null&&ge(_.ref,null,U,_,!0)},k=(_,C,E,B)=>{if(_==null)e(C.el=p(C.children),E,B);else{const H=C.el=_.el;C.children!==_.children&&r(H,C.children)}},j=(_,C,E,B)=>{_==null?e(C.el=c(C.children||""),E,B):C.el=_.el},v=(_,C,E,B)=>{[_.el,_.anchor]=f(_.children,C,E,B,_.el,_.anchor)},y=({el:_,anchor:C},E,B)=>{let H;for(;_&&_!==C;)H=h(_),e(_,E,B),_=H;e(C,E,B)},w=({el:_,anchor:C})=>{let E;for(;_&&_!==C;)E=h(_),t(_),_=E;t(C)},R=(_,C,E,B,H,U,g,b,P)=>{if(C.type==="svg"?g="svg":C.type==="math"&&(g="mathml"),_==null)N(C,E,B,H,U,g,b,P);else{const M=_.el&&_.el._isVueCE?_.el:null;try{M&&M._beginPatch(),T(_,C,H,U,g,b,P)}finally{M&&M._endPatch()}}},N=(_,C,E,B,H,U,g,b)=>{let P,M;const{props:Q,shapeFlag:G,transition:Z,dirs:ss}=_;if(P=_.el=o(_.type,U,Q&&Q.is,Q),G&8?i(P,_.children):G&16&&A(_.children,P,null,B,H,tl(_,U),g,b),ss&&Sa(_,null,B,"created"),S(P,_,_.scopeId,g,B),Q){for(const q in Q)q!=="value"&&!de(q)&&l(P,q,null,Q[q],U,B);"value"in Q&&l(P,"value",null,Q.value,U),(M=Q.onVnodeBeforeMount)&&In(M,B,_)}ss&&Sa(_,null,B,"beforeMount");const L=ld(H,Z);L&&Z.beforeEnter(P),e(P,C,E),((M=Q&&Q.onVnodeMounted)||L||ss)&&mn(()=>{M&&In(M,B,_),L&&Z.enter(P),ss&&Sa(_,null,B,"mounted")},H)},S=(_,C,E,B,H)=>{if(E&&d(_,E),B)for(let U=0;U<B.length;U++)d(_,B[U]);if(H){let U=H.subTree;if(C===U||ai(U.type)&&(U.ssContent===C||U.ssFallback===C)){const g=H.vnode;S(_,g,g.scopeId,g.slotScopeIds,H.parent)}}},A=(_,C,E,B,H,U,g,b,P=0)=>{for(let M=P;M<_.length;M++){const Q=_[M]=b?ja(_[M]):$n(_[M]);m(null,Q,C,E,B,H,U,g,b)}},T=(_,C,E,B,H,U,g)=>{const b=C.el=_.el;let{patchFlag:P,dynamicChildren:M,dirs:Q}=C;P|=_.patchFlag&16;const G=_.props||As,Z=C.props||As;let ss;if(E&&Ea(E,!1),(ss=Z.onVnodeBeforeUpdate)&&In(ss,E,C,_),Q&&Sa(C,_,E,"beforeUpdate"),E&&Ea(E,!0),(G.innerHTML&&Z.innerHTML==null||G.textContent&&Z.textContent==null)&&i(b,""),M?z(_.dynamicChildren,M,b,E,B,tl(C,H),U):g||cs(_,C,b,null,E,B,tl(C,H),U,!1),P>0){if(P&16)I(b,G,Z,E,H);else if(P&2&&G.class!==Z.class&&l(b,"class",null,Z.class,H),P&4&&l(b,"style",G.style,Z.style,H),P&8){const L=C.dynamicProps;for(let q=0;q<L.length;q++){const ts=L[q],bs=G[ts],Vs=Z[ts];(Vs!==bs||ts==="value")&&l(b,ts,bs,Vs,H,E)}}P&1&&_.children!==C.children&&i(b,C.children)}else!g&&M==null&&I(b,G,Z,E,H);((ss=Z.onVnodeUpdated)||Q)&&mn(()=>{ss&&In(ss,E,C,_),Q&&Sa(C,_,E,"updated")},B)},z=(_,C,E,B,H,U,g)=>{for(let b=0;b<C.length;b++){const P=_[b],M=C[b],Q=P.el&&(P.type===js||!La(P,M)||P.shapeFlag&198)?u(P.el):E;m(P,M,Q,null,B,H,U,g,!0)}},I=(_,C,E,B,H)=>{if(C!==E){if(C!==As)for(const U in C)!de(U)&&!(U in E)&&l(_,U,C[U],null,H,B);for(const U in E){if(de(U))continue;const g=E[U],b=C[U];g!==b&&U!=="value"&&l(_,U,b,g,H,B)}"value"in E&&l(_,"value",C.value,E.value,H)}},D=(_,C,E,B,H,U,g,b,P)=>{const M=C.el=_?_.el:p(""),Q=C.anchor=_?_.anchor:p("");let{patchFlag:G,dynamicChildren:Z,slotScopeIds:ss}=C;ss&&(b=b?b.concat(ss):ss),_==null?(e(M,E,B),e(Q,E,B),A(C.children||[],E,Q,H,U,g,b,P)):G>0&&G&64&&Z&&_.dynamicChildren?(z(_.dynamicChildren,Z,E,H,U,g,b),(C.key!=null||H&&C===H.subTree)&&Eo(_,C,!0)):cs(_,C,E,Q,H,U,g,b,P)},W=(_,C,E,B,H,U,g,b,P)=>{C.slotScopeIds=b,_==null?C.shapeFlag&512?H.ctx.activate(C,E,B,g,P):ls(C,E,B,H,U,g,P):ms(_,C,P)},ls=(_,C,E,B,H,U,g)=>{const b=_.component=xd(_,B,H);if(ze(_)&&(b.ctx.renderer=as),Pd(b,!1,g),b.asyncDep){if(H&&H.registerDep(b,us,g),!_.el){const P=b.subTree=fs(cn);j(null,P,C,E),_.placeholder=P.el}}else us(b,_,C,E,H,U,g)},ms=(_,C,E)=>{const B=C.component=_.component;if(fd(_,C,E))if(B.asyncDep&&!B.asyncResolved){es(B,C,E);return}else B.next=C,B.update();else C.el=_.el,B.vnode=C},us=(_,C,E,B,H,U,g)=>{const b=()=>{if(_.isMounted){let{next:G,bu:Z,u:ss,parent:L,vnode:q}=_;{const Js=Zr(_);if(Js){G&&(G.el=q.el,es(_,G,g)),Js.asyncDep.then(()=>{_.isUnmounted||b()});return}}let ts=G,bs;Ea(_,!1),G?(G.el=q.el,es(_,G,g)):G=q,Z&&ot(Z),(bs=G.props&&G.props.onVnodeBeforeUpdate)&&In(bs,L,G,q),Ea(_,!0);const Vs=hp(_),un=_.subTree;_.subTree=Vs,m(un,Vs,u(un.el),F(un),_,H,U),G.el=Vs.el,ts===null&&gd(_,Vs.el),ss&&mn(ss,H),(bs=G.props&&G.props.onVnodeUpdated)&&mn(()=>In(bs,L,G,q),H)}else{let G;const{el:Z,props:ss}=C,{bm:L,m:q,parent:ts,root:bs,type:Vs}=_,un=Ya(C);Ea(_,!1),L&&ot(L),!un&&(G=ss&&ss.onVnodeBeforeMount)&&In(G,ts,C),Ea(_,!0);{bs.ce&&bs.ce._def.shadowRoot!==!1&&bs.ce._injectChildStyle(Vs);const Js=_.subTree=hp(_);m(null,Js,E,B,_,H,U),C.el=Js.el}if(q&&mn(q,H),!un&&(G=ss&&ss.onVnodeMounted)){const Js=C;mn(()=>In(G,ts,Js),H)}(C.shapeFlag&256||ts&&Ya(ts.vnode)&&ts.vnode.shapeFlag&256)&&_.a&&mn(_.a,H),_.isMounted=!0,C=E=B=null}};_.scope.on();const P=_.effect=new tr(b);_.scope.off();const M=_.update=P.run.bind(P),Q=_.job=P.runIfDirty.bind(P);Q.i=_,Q.id=_.uid,P.scheduler=()=>vo(Q),Ea(_,!0),M()},es=(_,C,E)=>{C.component=_;const B=_.vnode.props;_.vnode=C,_.next=null,Jh(_,C.props,B,E),ad(_,C.children,E),sa(),Jo(_),na()},cs=(_,C,E,B,H,U,g,b,P=!1)=>{const M=_&&_.children,Q=_?_.shapeFlag:0,G=C.children,{patchFlag:Z,shapeFlag:ss}=C;if(Z>0){if(Z&128){Bs(M,G,E,B,H,U,g,b,P);return}else if(Z&256){$s(M,G,E,B,H,U,g,b,P);return}}ss&8?(Q&16&&Ts(M,H,U),G!==M&&i(E,G)):Q&16?ss&16?Bs(M,G,E,B,H,U,g,b,P):Ts(M,H,U,!0):(Q&8&&i(E,""),ss&16&&A(G,E,B,H,U,g,b,P))},$s=(_,C,E,B,H,U,g,b,P)=>{_=_||Wa,C=C||Wa;const M=_.length,Q=C.length,G=Math.min(M,Q);let Z;for(Z=0;Z<G;Z++){const ss=C[Z]=P?ja(C[Z]):$n(C[Z]);m(_[Z],ss,E,null,H,U,g,b,P)}M>Q?Ts(_,H,U,!0,!1,G):A(C,E,B,H,U,g,b,P,G)},Bs=(_,C,E,B,H,U,g,b,P)=>{let M=0;const Q=C.length;let G=_.length-1,Z=Q-1;for(;M<=G&&M<=Z;){const ss=_[M],L=C[M]=P?ja(C[M]):$n(C[M]);if(La(ss,L))m(ss,L,E,null,H,U,g,b,P);else break;M++}for(;M<=G&&M<=Z;){const ss=_[G],L=C[Z]=P?ja(C[Z]):$n(C[Z]);if(La(ss,L))m(ss,L,E,null,H,U,g,b,P);else break;G--,Z--}if(M>G){if(M<=Z){const ss=Z+1,L=ss<Q?C[ss].el:B;for(;M<=Z;)m(null,C[M]=P?ja(C[M]):$n(C[M]),E,L,H,U,g,b,P),M++}}else if(M>Z)for(;M<=G;)os(_[M],H,U,!0),M++;else{const ss=M,L=M,q=new Map;for(M=L;M<=Z;M++){const tn=C[M]=P?ja(C[M]):$n(C[M]);tn.key!=null&&q.set(tn.key,M)}let ts,bs=0;const Vs=Z-L+1;let un=!1,Js=0;const yn=new Array(Vs);for(M=0;M<Vs;M++)yn[M]=0;for(M=ss;M<=G;M++){const tn=_[M];if(bs>=Vs){os(tn,H,U,!0);continue}let On;if(tn.key!=null)On=q.get(tn.key);else for(ts=L;ts<=Z;ts++)if(yn[ts-L]===0&&La(tn,C[ts])){On=ts;break}On===void 0?os(tn,H,U,!0):(yn[On-L]=M+1,On>=Js?Js=On:un=!0,m(tn,C[On],E,null,H,U,g,b,P),bs++)}const Ke=un?od(yn):Wa;for(ts=Ke.length-1,M=Vs-1;M>=0;M--){const tn=L+M,On=C[tn],Go=C[tn+1],Wo=tn+1<Q?Go.el||Go.placeholder:B;yn[M]===0?m(null,On,E,Wo,H,U,g,b,P):un&&(ts<0||M!==Ke[ts]?qs(On,E,Wo,2):ts--)}}},qs=(_,C,E,B,H=null)=>{const{el:U,type:g,transition:b,children:P,shapeFlag:M}=_;if(M&6){qs(_.component.subTree,C,E,B);return}if(M&128){_.suspense.move(C,E,B);return}if(M&64){g.move(_,C,E,as);return}if(g===js){e(U,C,E);for(let G=0;G<P.length;G++)qs(P[G],C,E,B);e(_.anchor,C,E);return}if(g===ll){y(_,C,E);return}if(B!==2&&M&1&&b)if(B===0)b.beforeEnter(U),e(U,C,E),mn(()=>b.enter(U),H);else{const{leave:G,delayLeave:Z,afterLeave:ss}=b,L=()=>{_.ctx.isUnmounted?t(U):e(U,C,E)},q=()=>{U._isLeaving&&U[Xn](!0),G(U,()=>{L(),ss&&ss()})};Z?Z(U,L,q):q()}else e(U,C,E)},os=(_,C,E,B=!1,H=!1)=>{const{type:U,props:g,ref:b,children:P,dynamicChildren:M,shapeFlag:Q,patchFlag:G,dirs:Z,cacheIndex:ss}=_;if(G===-2&&(H=!1),b!=null&&(sa(),ge(b,null,E,_,!0),na()),ss!=null&&(C.renderCache[ss]=void 0),Q&256){C.ctx.deactivate(_);return}const L=Q&1&&Z,q=!Ya(_);let ts;if(q&&(ts=g&&g.onVnodeBeforeUnmount)&&In(ts,C,_),Q&6)Ds(_.component,E,B);else{if(Q&128){_.suspense.unmount(E,B);return}L&&Sa(_,null,C,"beforeUnmount"),Q&64?_.type.remove(_,C,E,as,B):M&&!M.hasOnce&&(U!==js||G>0&&G&64)?Ts(M,C,E,!1,!0):(U===js&&G&384||!H&&Q&16)&&Ts(P,C,E),B&&gs(_)}(q&&(ts=g&&g.onVnodeUnmounted)||L)&&mn(()=>{ts&&In(ts,C,_),L&&Sa(_,null,C,"unmounted")},E)},gs=_=>{const{type:C,el:E,anchor:B,transition:H}=_;if(C===js){_s(E,B);return}if(C===ll){w(_);return}const U=()=>{t(E),H&&!H.persisted&&H.afterLeave&&H.afterLeave()};if(_.shapeFlag&1&&H&&!H.persisted){const{leave:g,delayLeave:b}=H,P=()=>g(E,U);b?b(_.el,U,P):P()}else U()},_s=(_,C)=>{let E;for(;_!==C;)E=h(_),t(_),_=E;t(C)},Ds=(_,C,E)=>{const{bum:B,scope:H,job:U,subTree:g,um:b,m:P,a:M}=_;up(P),up(M),B&&ot(B),H.stop(),U&&(U.flags|=8,os(g,_,C,E)),b&&mn(b,C),mn(()=>{_.isUnmounted=!0},C)},Ts=(_,C,E,B=!1,H=!1,U=0)=>{for(let g=U;g<_.length;g++)os(_[g],C,E,B,H)},F=_=>{if(_.shapeFlag&6)return F(_.component.subTree);if(_.shapeFlag&128)return _.suspense.next();const C=h(_.anchor||_.el),E=C&&C[Pr];return E?h(E):C};let X=!1;const K=(_,C,E)=>{_==null?C._vnode&&os(C._vnode,null,null,!0):m(C._vnode||null,_,C,null,null,null,E),C._vnode=_,X||(X=!0,Jo(),Cr(),X=!1)},as={p:m,um:os,m:qs,r:gs,mt:ls,mc:A,pc:cs,pbc:z,n:F,o:s};return{render:K,hydrate:void 0,createApp:Yh(K)}}function tl({type:s,props:n},a){return a==="svg"&&s==="foreignObject"||a==="mathml"&&s==="annotation-xml"&&n&&n.encoding&&n.encoding.includes("html")?void 0:a}function Ea({effect:s,job:n},a){a?(s.flags|=32,n.flags|=4):(s.flags&=-33,n.flags&=-5)}function ld(s,n){return(!s||s&&!s.pendingBranch)&&n&&!n.persisted}function Eo(s,n,a=!1){const e=s.children,t=n.children;if(rs(e)&&rs(t))for(let l=0;l<e.length;l++){const o=e[l];let p=t[l];p.shapeFlag&1&&!p.dynamicChildren&&((p.patchFlag<=0||p.patchFlag===32)&&(p=t[l]=ja(t[l]),p.el=o.el),!a&&p.patchFlag!==-2&&Eo(o,p)),p.type===Ue&&p.patchFlag!==-1&&(p.el=o.el),p.type===cn&&!p.el&&(p.el=o.el)}}function od(s){const n=s.slice(),a=[0];let e,t,l,o,p;const c=s.length;for(e=0;e<c;e++){const r=s[e];if(r!==0){if(t=a[a.length-1],s[t]<r){n[e]=t,a.push(e);continue}for(l=0,o=a.length-1;l<o;)p=l+o>>1,s[a[p]]<r?l=p+1:o=p;r<s[a[l]]&&(l>0&&(n[e]=a[l-1]),a[l]=e)}}for(l=a.length,o=a[l-1];l-- >0;)a[l]=o,o=n[o];return a}function Zr(s){const n=s.subTree.component;if(n)return n.asyncDep&&!n.asyncResolved?n:Zr(n)}function up(s){if(s)for(let n=0;n<s.length;n++)s[n].flags|=8}const pd=Symbol.for("v-scx"),cd=()=>_n(pd);function rd(s,n){return To(s,null,n)}function Ws(s,n,a){return To(s,n,a)}function To(s,n,a=As){const{immediate:e,deep:t,flush:l,once:o}=a,p=Ys({},a),c=n&&e||!n&&l!=="post";let r;if(Ja){if(l==="sync"){const d=cd();r=d.__watcherHandles||(d.__watcherHandles=[])}else if(!c){const d=()=>{};return d.stop=qn,d.resume=qn,d.pause=qn,d}}const i=sn;p.call=(d,f,m)=>Ln(d,i,f,m);let u=!1;l==="post"?p.scheduler=d=>{mn(d,i&&i.suspense)}:l!=="sync"&&(u=!0,p.scheduler=(d,f)=>{f?d():vo(d)}),p.augmentJob=d=>{n&&(d.flags|=4),u&&(d.flags|=2,i&&(d.id=i.uid,d.i=i))};const h=Ch(s,n,p);return Ja&&(r?r.push(h):c&&h()),h}function id(s,n,a){const e=this.proxy,t=Us(s)?s.includes(".")?si(e,s):()=>e[s]:s.bind(e,e);let l;ds(n)?l=n:(l=n.handler,a=n);const o=He(this),p=To(t,l.bind(e),a);return o(),p}function si(s,n){const a=n.split(".");return()=>{let e=s;for(let t=0;t<a.length&&e;t++)e=e[a[t]];return e}}const ud=(s,n)=>n==="modelValue"||n==="model-value"?s.modelModifiers:s[`${n}Modifiers`]||s[`${Tn(n)}Modifiers`]||s[`${Ca(n)}Modifiers`];function hd(s,n,...a){if(s.isUnmounted)return;const e=s.vnode.props||As;let t=a;const l=n.startsWith("update:"),o=l&&ud(e,n.slice(7));o&&(o.trim&&(t=a.map(i=>Us(i)?i.trim():i)),o.number&&(t=a.map(co)));let p,c=e[p=Qt(n)]||e[p=Qt(Tn(n))];!c&&l&&(c=e[p=Qt(Ca(n))]),c&&Ln(c,s,6,t);const r=e[p+"Once"];if(r){if(!s.emitted)s.emitted={};else if(s.emitted[p])return;s.emitted[p]=!0,Ln(r,s,6,t)}}const dd=new WeakMap;function ni(s,n,a=!1){const e=a?dd:n.emitsCache,t=e.get(s);if(t!==void 0)return t;const l=s.emits;let o={},p=!1;if(!ds(s)){const c=r=>{const i=ni(r,n,!0);i&&(p=!0,Ys(o,i))};!a&&n.mixins.length&&n.mixins.forEach(c),s.extends&&c(s.extends),s.mixins&&s.mixins.forEach(c)}return!l&&!p?(Ms(s)&&e.set(s,null),null):(rs(l)?l.forEach(c=>o[c]=null):Ys(o,l),Ms(s)&&e.set(s,o),o)}function Ft(s,n){return!s||!Rt(n)?!1:(n=n.slice(2).replace(/Once$/,""),Ss(s,n[0].toLowerCase()+n.slice(1))||Ss(s,Ca(n))||Ss(s,n))}function hp(s){const{type:n,vnode:a,proxy:e,withProxy:t,propsOptions:[l],slots:o,attrs:p,emit:c,render:r,renderCache:i,props:u,data:h,setupState:d,ctx:f,inheritAttrs:m}=s,k=bt(s);let j,v;try{if(a.shapeFlag&4){const w=t||e,R=w;j=$n(r.call(R,w,i,u,d,h,f)),v=p}else{const w=n;j=$n(w.length>1?w(u,{attrs:p,slots:o,emit:c}):w(u,null)),v=n.props?p:md(p)}}catch(w){_e.length=0,qe(w,s,1),j=fs(cn)}let y=j;if(v&&m!==!1){const w=Object.keys(v),{shapeFlag:R}=y;w.length&&R&7&&(l&&w.some(lo)&&(v=jd(v,l)),y=va(y,v,!1,!0))}return a.dirs&&(y=va(y,null,!1,!0),y.dirs=y.dirs?y.dirs.concat(a.dirs):a.dirs),a.transition&&Se(y,a.transition),j=y,bt(k),j}const md=s=>{let n;for(const a in s)(a==="class"||a==="style"||Rt(a))&&((n||(n={}))[a]=s[a]);return n},jd=(s,n)=>{const a={};for(const e in s)(!lo(e)||!(e.slice(9)in n))&&(a[e]=s[e]);return a};function fd(s,n,a){const{props:e,children:t,component:l}=s,{props:o,children:p,patchFlag:c}=n,r=l.emitsOptions;if(n.dirs||n.transition)return!0;if(a&&c>=0){if(c&1024)return!0;if(c&16)return e?dp(e,o,r):!!o;if(c&8){const i=n.dynamicProps;for(let u=0;u<i.length;u++){const h=i[u];if(o[h]!==e[h]&&!Ft(r,h))return!0}}}else return(t||p)&&(!p||!p.$stable)?!0:e===o?!1:e?o?dp(e,o,r):!0:!!o;return!1}function dp(s,n,a){const e=Object.keys(n);if(e.length!==Object.keys(s).length)return!0;for(let t=0;t<e.length;t++){const l=e[t];if(n[l]!==s[l]&&!Ft(a,l))return!0}return!1}function gd({vnode:s,parent:n},a){for(;n;){const e=n.subTree;if(e.suspense&&e.suspense.activeBranch===s&&(e.el=s.el),e===s)(s=n.vnode).el=a,n=n.parent;else break}}const ai=s=>s.__isSuspense;function bd(s,n){n&&n.pendingBranch?rs(s)?n.effects.push(...s):n.effects.push(s):Ph(s)}const js=Symbol.for("v-fgt"),Ue=Symbol.for("v-txt"),cn=Symbol.for("v-cmt"),ll=Symbol.for("v-stc"),_e=[];let wn=null;function O(s=!1){_e.push(wn=s?null:[])}function _d(){_e.pop(),wn=_e[_e.length-1]||null}let Ee=1;function vt(s,n=!1){Ee+=s,s<0&&wn&&n&&(wn.hasOnce=!0)}function ei(s){return s.dynamicChildren=Ee>0?wn||Wa:null,_d(),Ee>0&&wn&&wn.push(s),s}function $(s,n,a,e,t,l){return ei(x(s,n,a,e,t,l,!0))}function rn(s,n,a,e,t){return ei(fs(s,n,a,e,t,!0))}function Te(s){return s?s.__v_isVNode===!0:!1}function La(s,n){return s.type===n.type&&s.key===n.key}const ti=({key:s})=>s??null,it=({ref:s,ref_key:n,ref_for:a})=>(typeof s=="number"&&(s=""+s),s!=null?Us(s)||zs(s)||ds(s)?{i:nn,r:s,k:n,f:!!a}:s:null);function x(s,n=null,a=null,e=0,t=null,l=s===js?0:1,o=!1,p=!1){const c={__v_isVNode:!0,__v_skip:!0,type:s,props:n,key:n&&ti(n),ref:n&&it(n),scopeId:xr,slotScopeIds:null,children:a,component:null,suspense:null,ssContent:null,ssFallback:null,dirs:null,transition:null,el:null,anchor:null,target:null,targetStart:null,targetAnchor:null,staticCount:0,shapeFlag:l,patchFlag:e,dynamicProps:t,dynamicChildren:null,appContext:null,ctx:nn};return p?(Ao(c,a),l&128&&s.normalize(c)):a&&(c.shapeFlag|=Us(a)?8:16),Ee>0&&!o&&wn&&(c.patchFlag>0||l&6)&&c.patchFlag!==32&&wn.push(c),c}const fs=yd;function yd(s,n=null,a=null,e=0,t=null,l=!1){if((!s||s===Fr)&&(s=cn),Te(s)){const p=va(s,n,!0);return a&&Ao(p,a),Ee>0&&!l&&wn&&(p.shapeFlag&6?wn[wn.indexOf(s)]=p:wn.push(p)),p.patchFlag=-2,p}if(Rd(s)&&(s=s.__vccOpts),n){n=vd(n);let{class:p,style:c}=n;p&&!Us(p)&&(n.class=Fs(p)),Ms(c)&&(bo(c)&&!rs(c)&&(c=Ys({},c)),n.style=oa(c))}const o=Us(s)?1:ai(s)?128:Sr(s)?64:Ms(s)?4:ds(s)?2:0;return x(s,n,a,e,t,o,l,!0)}function vd(s){return s?bo(s)||Wr(s)?Ys({},s):s:null}function va(s,n,a=!1,e=!1){const{props:t,ref:l,patchFlag:o,children:p,transition:c}=s,r=n?wd(t||{},n):t,i={__v_isVNode:!0,__v_skip:!0,type:s.type,props:r,key:r&&ti(r),ref:n&&n.ref?a&&l?rs(l)?l.concat(it(n)):[l,it(n)]:it(n):l,scopeId:s.scopeId,slotScopeIds:s.slotScopeIds,children:p,target:s.target,targetStart:s.targetStart,targetAnchor:s.targetAnchor,staticCount:s.staticCount,shapeFlag:s.shapeFlag,patchFlag:n&&s.type!==js?o===-1?16:o|16:o,dynamicProps:s.dynamicProps,dynamicChildren:s.dynamicChildren,appContext:s.appContext,dirs:s.dirs,transition:c,component:s.component,suspense:s.suspense,ssContent:s.ssContent&&va(s.ssContent),ssFallback:s.ssFallback&&va(s.ssFallback),placeholder:s.placeholder,el:s.el,anchor:s.anchor,ctx:s.ctx,ce:s.ce};return c&&e&&Se(i,c.clone(i)),i}function ks(s=" ",n=0){return fs(Ue,null,s,n)}function is(s="",n=!1){return n?(O(),rn(cn,null,s)):fs(cn,null,s)}function $n(s){return s==null||typeof s=="boolean"?fs(cn):rs(s)?fs(js,null,s.slice()):Te(s)?ja(s):fs(Ue,null,String(s))}function ja(s){return s.el===null&&s.patchFlag!==-1||s.memo?s:va(s)}function Ao(s,n){let a=0;const{shapeFlag:e}=s;if(n==null)n=null;else if(rs(n))a=16;else if(typeof n=="object")if(e&65){const t=n.default;t&&(t._c&&(t._d=!1),Ao(s,t()),t._c&&(t._d=!0));return}else{a=32;const t=n._;!t&&!Wr(n)?n._ctx=nn:t===3&&nn&&(nn.slots._===1?n._=1:(n._=2,s.patchFlag|=1024))}else ds(n)?(n={default:n,_ctx:nn},a=32):(n=String(n),e&64?(a=16,n=[ks(n)]):a=8);s.children=n,s.shapeFlag|=a}function wd(...s){const n={};for(let a=0;a<s.length;a++){const e=s[a];for(const t in e)if(t==="class")n.class!==e.class&&(n.class=Fs([n.class,e.class]));else if(t==="style")n.style=oa([n.style,e.style]);else if(Rt(t)){const l=n[t],o=e[t];o&&l!==o&&!(rs(l)&&l.includes(o))&&(n[t]=l?[].concat(l,o):o)}else t!==""&&(n[t]=e[t])}return n}function In(s,n,a,e=null){Ln(s,n,7,[a,e])}const Cd=Ur();let kd=0;function xd(s,n,a){const e=s.type,t=(n?n.appContext:s.appContext)||Cd,l={uid:kd++,vnode:s,type:e,parent:n,appContext:t,root:null,next:null,subTree:null,effect:null,update:null,job:null,scope:new ar(!0),render:null,proxy:null,exposed:null,exposeProxy:null,withProxy:null,provides:n?n.provides:Object.create(t.provides),ids:n?n.ids:["",0,0],accessCache:null,renderCache:[],components:null,directives:null,propsOptions:Xr(e,t),emitsOptions:ni(e,t),emit:null,emitted:null,propsDefaults:As,inheritAttrs:e.inheritAttrs,ctx:As,data:As,props:As,attrs:As,slots:As,refs:As,setupState:As,setupContext:null,suspense:a,suspenseId:a?a.pendingId:0,asyncDep:null,asyncResolved:!1,isMounted:!1,isUnmounted:!1,isDeactivated:!1,bc:null,c:null,bm:null,m:null,bu:null,u:null,um:null,bum:null,da:null,a:null,rtg:null,rtc:null,ec:null,sp:null};return l.ctx={_:l},l.root=n?n.root:l,l.emit=hd.bind(null,l),s.ce&&s.ce(l),l}let sn=null;const ca=()=>sn||nn;let wt,Ll;{const s=Ot(),n=(a,e)=>{let t;return(t=s[a])||(t=s[a]=[]),t.push(e),l=>{t.length>1?t.forEach(o=>o(l)):t[0](l)}};wt=n("__VUE_INSTANCE_SETTERS__",a=>sn=a),Ll=n("__VUE_SSR_SETTERS__",a=>Ja=a)}const He=s=>{const n=sn;return wt(s),s.scope.on(),()=>{s.scope.off(),wt(n)}},mp=()=>{sn&&sn.scope.off(),wt(null)};function li(s){return s.vnode.shapeFlag&4}let Ja=!1;function Pd(s,n=!1,a=!1){n&&Ll(n);const{props:e,children:t}=s.vnode,l=li(s);Qh(s,e,l,n),nd(s,t,a||n);const o=l?Sd(s,n):void 0;return n&&Ll(!1),o}function Sd(s,n){const a=s.type;s.accessCache=Object.create(null),s.proxy=new Proxy(s.ctx,Uh);const{setup:e}=a;if(e){sa();const t=s.setupContext=e.length>1?Td(s):null,l=He(s),o=Be(e,s,0,[s.props,t]),p=Xc(o);if(na(),l(),(p||s.sp)&&!Ya(s)&&wo(s),p){if(o.then(mp,mp),n)return o.then(c=>{jp(s,c)}).catch(c=>{qe(c,s,0)});s.asyncDep=o}else jp(s,o)}else oi(s)}function jp(s,n,a){ds(n)?s.type.__ssrInlineRender?s.ssrRender=n:s.render=n:Ms(n)&&(s.setupState=yr(n)),oi(s)}function oi(s,n,a){const e=s.type;s.render||(s.render=e.render||qn);{const t=He(s);sa();try{Hh(s)}finally{na(),t()}}}const Ed={get(s,n){return pn(s,"get",""),s[n]}};function Td(s){const n=a=>{s.exposed=a||{}};return{attrs:new Proxy(s.attrs,Ed),slots:s.slots,emit:s.emit,expose:n}}function $t(s){return s.exposed?s.exposeProxy||(s.exposeProxy=new Proxy(yr(_o(s.exposed)),{get(n,a){if(a in n)return n[a];if(a in be)return be[a](s)},has(n,a){return a in n||a in be}})):s.proxy}function Ad(s,n=!0){return ds(s)?s.displayName||s.name:s.name||n&&s.__name}function Rd(s){return ds(s)&&"__vccOpts"in s}const J=(s,n)=>vh(s,n,Ja);function Ve(s,n,a){try{vt(-1);const e=arguments.length;return e===2?Ms(n)&&!rs(n)?Te(n)?fs(s,null,[n]):fs(s,n):fs(s,null,n):(e>3?a=Array.prototype.slice.call(arguments,2):e===3&&Te(a)&&(a=[a]),fs(s,n,a))}finally{vt(1)}}const Ld="3.5.24";/**
* @vue/runtime-dom v3.5.24
* (c) 2018-present Yuxi (Evan) You and Vue contributors
* @license MIT
**/let Ml;const fp=typeof window<"u"&&window.trustedTypes;if(fp)try{Ml=fp.createPolicy("vue",{createHTML:s=>s})}catch{}const pi=Ml?s=>Ml.createHTML(s):s=>s,Md="http://www.w3.org/2000/svg",Dd="http://www.w3.org/1998/Math/MathML",Kn=typeof document<"u"?document:null,gp=Kn&&Kn.createElement("template"),Od={insert:(s,n,a)=>{n.insertBefore(s,a||null)},remove:s=>{const n=s.parentNode;n&&n.removeChild(s)},createElement:(s,n,a,e)=>{const t=n==="svg"?Kn.createElementNS(Md,s):n==="mathml"?Kn.createElementNS(Dd,s):a?Kn.createElement(s,{is:a}):Kn.createElement(s);return s==="select"&&e&&e.multiple!=null&&t.setAttribute("multiple",e.multiple),t},createText:s=>Kn.createTextNode(s),createComment:s=>Kn.createComment(s),setText:(s,n)=>{s.nodeValue=n},setElementText:(s,n)=>{s.textContent=n},parentNode:s=>s.parentNode,nextSibling:s=>s.nextSibling,querySelector:s=>Kn.querySelector(s),setScopeId(s,n){s.setAttribute(n,"")},insertStaticContent(s,n,a,e,t,l){const o=a?a.previousSibling:n.lastChild;if(t&&(t===l||t.nextSibling))for(;n.insertBefore(t.cloneNode(!0),a),!(t===l||!(t=t.nextSibling)););else{gp.innerHTML=pi(e==="svg"?`<svg>${s}</svg>`:e==="mathml"?`<math>${s}</math>`:s);const p=gp.content;if(e==="svg"||e==="mathml"){const c=p.firstChild;for(;c.firstChild;)p.appendChild(c.firstChild);p.removeChild(c)}n.insertBefore(p,a)}return[o?o.nextSibling:n.firstChild,a?a.previousSibling:n.lastChild]}},ia="transition",pe="animation",Ae=Symbol("_vtc"),ci={name:String,type:String,css:{type:Boolean,default:!0},duration:[String,Number,Object],enterFromClass:String,enterActiveClass:String,enterToClass:String,appearFromClass:String,appearActiveClass:String,appearToClass:String,leaveFromClass:String,leaveActiveClass:String,leaveToClass:String},Id=Ys({},Ar,ci),Nd=s=>(s.displayName="Transition",s.props=Id,s),Fd=Nd((s,{slots:n})=>Ve(Ah,$d(s),n)),Ta=(s,n=[])=>{rs(s)?s.forEach(a=>a(...n)):s&&s(...n)},bp=s=>s?rs(s)?s.some(n=>n.length>1):s.length>1:!1;function $d(s){const n={};for(const D in s)D in ci||(n[D]=s[D]);if(s.css===!1)return n;const{name:a="v",type:e,duration:t,enterFromClass:l=`${a}-enter-from`,enterActiveClass:o=`${a}-enter-active`,enterToClass:p=`${a}-enter-to`,appearFromClass:c=l,appearActiveClass:r=o,appearToClass:i=p,leaveFromClass:u=`${a}-leave-from`,leaveActiveClass:h=`${a}-leave-active`,leaveToClass:d=`${a}-leave-to`}=s,f=Bd(t),m=f&&f[0],k=f&&f[1],{onBeforeEnter:j,onEnter:v,onEnterCancelled:y,onLeave:w,onLeaveCancelled:R,onBeforeAppear:N=j,onAppear:S=v,onAppearCancelled:A=y}=n,T=(D,W,ls,ms)=>{D._enterCancelled=ms,Aa(D,W?i:p),Aa(D,W?r:o),ls&&ls()},z=(D,W)=>{D._isLeaving=!1,Aa(D,u),Aa(D,d),Aa(D,h),W&&W()},I=D=>(W,ls)=>{const ms=D?S:v,us=()=>T(W,D,ls);Ta(ms,[W,us]),_p(()=>{Aa(W,D?c:l),Vn(W,D?i:p),bp(ms)||yp(W,e,m,us)})};return Ys(n,{onBeforeEnter(D){Ta(j,[D]),Vn(D,l),Vn(D,o)},onBeforeAppear(D){Ta(N,[D]),Vn(D,c),Vn(D,r)},onEnter:I(!1),onAppear:I(!0),onLeave(D,W){D._isLeaving=!0;const ls=()=>z(D,W);Vn(D,u),D._enterCancelled?(Vn(D,h),Cp(D)):(Cp(D),Vn(D,h)),_p(()=>{D._isLeaving&&(Aa(D,u),Vn(D,d),bp(w)||yp(D,e,k,ls))}),Ta(w,[D,ls])},onEnterCancelled(D){T(D,!1,void 0,!0),Ta(y,[D])},onAppearCancelled(D){T(D,!0,void 0,!0),Ta(A,[D])},onLeaveCancelled(D){z(D),Ta(R,[D])}})}function Bd(s){if(s==null)return null;if(Ms(s))return[ol(s.enter),ol(s.leave)];{const n=ol(s);return[n,n]}}function ol(s){return zu(s)}function Vn(s,n){n.split(/\s+/).forEach(a=>a&&s.classList.add(a)),(s[Ae]||(s[Ae]=new Set)).add(n)}function Aa(s,n){n.split(/\s+/).forEach(e=>e&&s.classList.remove(e));const a=s[Ae];a&&(a.delete(n),a.size||(s[Ae]=void 0))}function _p(s){requestAnimationFrame(()=>{requestAnimationFrame(s)})}let qd=0;function yp(s,n,a,e){const t=s._endId=++qd,l=()=>{t===s._endId&&e()};if(a!=null)return setTimeout(l,a);const{type:o,timeout:p,propCount:c}=zd(s,n);if(!o)return e();const r=o+"end";let i=0;const u=()=>{s.removeEventListener(r,h),l()},h=d=>{d.target===s&&++i>=c&&u()};setTimeout(()=>{i<c&&u()},p+1),s.addEventListener(r,h)}function zd(s,n){const a=window.getComputedStyle(s),e=f=>(a[f]||"").split(", "),t=e(`${ia}Delay`),l=e(`${ia}Duration`),o=vp(t,l),p=e(`${pe}Delay`),c=e(`${pe}Duration`),r=vp(p,c);let i=null,u=0,h=0;n===ia?o>0&&(i=ia,u=o,h=l.length):n===pe?r>0&&(i=pe,u=r,h=c.length):(u=Math.max(o,r),i=u>0?o>r?ia:pe:null,h=i?i===ia?l.length:c.length:0);const d=i===ia&&/\b(?:transform|all)(?:,|$)/.test(e(`${ia}Property`).toString());return{type:i,timeout:u,propCount:h,hasTransform:d}}function vp(s,n){for(;s.length<n.length;)s=s.concat(s);return Math.max(...n.map((a,e)=>wp(a)+wp(s[e])))}function wp(s){return s==="auto"?0:Number(s.slice(0,-1).replace(",","."))*1e3}function Cp(s){return(s?s.ownerDocument:document).body.offsetHeight}function Ud(s,n,a){const e=s[Ae];e&&(n=(n?[n,...e]:[...e]).join(" ")),n==null?s.removeAttribute("class"):a?s.setAttribute("class",n):s.className=n}const Ct=Symbol("_vod"),ri=Symbol("_vsh"),ii={name:"show",beforeMount(s,{value:n},{transition:a}){s[Ct]=s.style.display==="none"?"":s.style.display,a&&n?a.beforeEnter(s):ce(s,n)},mounted(s,{value:n},{transition:a}){a&&n&&a.enter(s)},updated(s,{value:n,oldValue:a},{transition:e}){!n!=!a&&(e?n?(e.beforeEnter(s),ce(s,!0),e.enter(s)):e.leave(s,()=>{ce(s,!1)}):ce(s,n))},beforeUnmount(s,{value:n}){ce(s,n)}};function ce(s,n){s.style.display=n?s[Ct]:"none",s[ri]=!n}const Hd=Symbol(""),Vd=/(?:^|;)\s*display\s*:/;function Gd(s,n,a){const e=s.style,t=Us(a);let l=!1;if(a&&!t){if(n)if(Us(n))for(const o of n.split(";")){const p=o.slice(0,o.indexOf(":")).trim();a[p]==null&&ut(e,p,"")}else for(const o in n)a[o]==null&&ut(e,o,"");for(const o in a)o==="display"&&(l=!0),ut(e,o,a[o])}else if(t){if(n!==a){const o=e[Hd];o&&(a+=";"+o),e.cssText=a,l=Vd.test(a)}}else n&&s.removeAttribute("style");Ct in s&&(s[Ct]=l?e.display:"",s[ri]&&(e.display="none"))}const kp=/\s*!important$/;function ut(s,n,a){if(rs(a))a.forEach(e=>ut(s,n,e));else if(a==null&&(a=""),n.startsWith("--"))s.setProperty(n,a);else{const e=Wd(s,n);kp.test(a)?s.setProperty(Ca(e),a.replace(kp,""),"important"):s[e]=a}}const xp=["Webkit","Moz","ms"],pl={};function Wd(s,n){const a=pl[n];if(a)return a;let e=Tn(n);if(e!=="filter"&&e in s)return pl[n]=e;e=Dt(e);for(let t=0;t<xp.length;t++){const l=xp[t]+e;if(l in s)return pl[n]=l}return n}const Pp="http://www.w3.org/1999/xlink";function Sp(s,n,a,e,t,l=Ku(n)){e&&n.startsWith("xlink:")?a==null?s.removeAttributeNS(Pp,n.slice(6,n.length)):s.setAttributeNS(Pp,n,a):a==null||l&&!Zc(a)?s.removeAttribute(n):s.setAttribute(n,l?"":la(a)?String(a):a)}function Ep(s,n,a,e,t){if(n==="innerHTML"||n==="textContent"){a!=null&&(s[n]=n==="innerHTML"?pi(a):a);return}const l=s.tagName;if(n==="value"&&l!=="PROGRESS"&&!l.includes("-")){const p=l==="OPTION"?s.getAttribute("value")||"":s.value,c=a==null?s.type==="checkbox"?"on":"":String(a);(p!==c||!("_value"in s))&&(s.value=c),a==null&&s.removeAttribute(n),s._value=a;return}let o=!1;if(a===""||a==null){const p=typeof s[n];p==="boolean"?a=Zc(a):a==null&&p==="string"?(a="",o=!0):p==="number"&&(a=0,o=!0)}try{s[n]=a}catch{}o&&s.removeAttribute(t||n)}function Ha(s,n,a,e){s.addEventListener(n,a,e)}function Kd(s,n,a,e){s.removeEventListener(n,a,e)}const Tp=Symbol("_vei");function Xd(s,n,a,e,t=null){const l=s[Tp]||(s[Tp]={}),o=l[n];if(e&&o)o.value=e;else{const[p,c]=Yd(n);if(e){const r=l[n]=Zd(e,t);Ha(s,p,r,c)}else o&&(Kd(s,p,o,c),l[n]=void 0)}}const Ap=/(?:Once|Passive|Capture)$/;function Yd(s){let n;if(Ap.test(s)){n={};let e;for(;e=s.match(Ap);)s=s.slice(0,s.length-e[0].length),n[e[0].toLowerCase()]=!0}return[s[2]===":"?s.slice(3):Ca(s.slice(2)),n]}let cl=0;const Qd=Promise.resolve(),Jd=()=>cl||(Qd.then(()=>cl=0),cl=Date.now());function Zd(s,n){const a=e=>{if(!e._vts)e._vts=Date.now();else if(e._vts<=a.attached)return;Ln(sm(e,a.value),n,5,[e])};return a.value=s,a.attached=Jd(),a}function sm(s,n){if(rs(n)){const a=s.stopImmediatePropagation;return s.stopImmediatePropagation=()=>{a.call(s),s._stopped=!0},n.map(e=>t=>!t._stopped&&e&&e(t))}else return n}const Rp=s=>s.charCodeAt(0)===111&&s.charCodeAt(1)===110&&s.charCodeAt(2)>96&&s.charCodeAt(2)<123,nm=(s,n,a,e,t,l)=>{const o=t==="svg";n==="class"?Ud(s,e,o):n==="style"?Gd(s,a,e):Rt(n)?lo(n)||Xd(s,n,a,e,l):(n[0]==="."?(n=n.slice(1),!0):n[0]==="^"?(n=n.slice(1),!1):am(s,n,e,o))?(Ep(s,n,e),!s.tagName.includes("-")&&(n==="value"||n==="checked"||n==="selected")&&Sp(s,n,e,o,l,n!=="value")):s._isVueCE&&(/[A-Z]/.test(n)||!Us(e))?Ep(s,Tn(n),e,l,n):(n==="true-value"?s._trueValue=e:n==="false-value"&&(s._falseValue=e),Sp(s,n,e,o))};function am(s,n,a,e){if(e)return!!(n==="innerHTML"||n==="textContent"||n in s&&Rp(n)&&ds(a));if(n==="spellcheck"||n==="draggable"||n==="translate"||n==="autocorrect"||n==="sandbox"&&s.tagName==="IFRAME"||n==="form"||n==="list"&&s.tagName==="INPUT"||n==="type"&&s.tagName==="TEXTAREA")return!1;if(n==="width"||n==="height"){const t=s.tagName;if(t==="IMG"||t==="VIDEO"||t==="CANVAS"||t==="SOURCE")return!1}return Rp(n)&&Us(a)?!1:n in s}const Lp=s=>{const n=s.props["onUpdate:modelValue"]||!1;return rs(n)?a=>ot(n,a):n};function em(s){s.target.composing=!0}function Mp(s){const n=s.target;n.composing&&(n.composing=!1,n.dispatchEvent(new Event("input")))}const rl=Symbol("_assign");function Dp(s,n,a){return n&&(s=s.trim()),a&&(s=co(s)),s}const k2={created(s,{modifiers:{lazy:n,trim:a,number:e}},t){s[rl]=Lp(t);const l=e||t.props&&t.props.type==="number";Ha(s,n?"change":"input",o=>{o.target.composing||s[rl](Dp(s.value,a,l))}),(a||l)&&Ha(s,"change",()=>{s.value=Dp(s.value,a,l)}),n||(Ha(s,"compositionstart",em),Ha(s,"compositionend",Mp),Ha(s,"change",Mp))},mounted(s,{value:n}){s.value=n??""},beforeUpdate(s,{value:n,oldValue:a,modifiers:{lazy:e,trim:t,number:l}},o){if(s[rl]=Lp(o),s.composing)return;const p=(l||s.type==="number")&&!/^0\d/.test(s.value)?co(s.value):s.value,c=n??"";p!==c&&(document.activeElement===s&&s.type!=="range"&&(e&&n===a||t&&s.value.trim()===c)||(s.value=c))}},tm=["ctrl","shift","alt","meta"],lm={stop:s=>s.stopPropagation(),prevent:s=>s.preventDefault(),self:s=>s.target!==s.currentTarget,ctrl:s=>!s.ctrlKey,shift:s=>!s.shiftKey,alt:s=>!s.altKey,meta:s=>!s.metaKey,left:s=>"button"in s&&s.button!==0,middle:s=>"button"in s&&s.button!==1,right:s=>"button"in s&&s.button!==2,exact:(s,n)=>tm.some(a=>s[`${a}Key`]&&!n.includes(a))},Bn=(s,n)=>{const a=s._withMods||(s._withMods={}),e=n.join(".");return a[e]||(a[e]=((t,...l)=>{for(let o=0;o<n.length;o++){const p=lm[n[o]];if(p&&p(t,n))return}return s(t,...l)}))},om={esc:"escape",space:" ",up:"arrow-up",left:"arrow-left",right:"arrow-right",down:"arrow-down",delete:"backspace"},Op=(s,n)=>{const a=s._withKeys||(s._withKeys={}),e=n.join(".");return a[e]||(a[e]=(t=>{if(!("key"in t))return;const l=Ca(t.key);if(n.some(o=>o===l||om[o]===l))return s(t)}))},pm=Ys({patchProp:nm},Od);let Ip;function cm(){return Ip||(Ip=ed(pm))}const rm=((...s)=>{const n=cm().createApp(...s),{mount:a}=n;return n.mount=e=>{const t=um(e);if(!t)return;const l=n._component;!ds(l)&&!l.render&&!l.template&&(l.template=t.innerHTML),t.nodeType===1&&(t.textContent="");const o=a(t,!1,im(t));return t instanceof Element&&(t.removeAttribute("v-cloak"),t.setAttribute("data-v-app","")),o},n});function im(s){if(s instanceof SVGElement)return"svg";if(typeof MathMLElement=="function"&&s instanceof MathMLElement)return"mathml"}function um(s){return Us(s)?document.querySelector(s):s}/*!
 * pinia v3.0.3
 * (c) 2025 Eduardo San Martin Morote
 * @license MIT
 */let ui;const Bt=s=>ui=s,hi=Symbol();function Dl(s){return s&&typeof s=="object"&&Object.prototype.toString.call(s)==="[object Object]"&&typeof s.toJSON!="function"}var ye;(function(s){s.direct="direct",s.patchObject="patch object",s.patchFunction="patch function"})(ye||(ye={}));function hm(){const s=ro(!0),n=s.run(()=>ps({}));let a=[],e=[];const t=_o({install(l){Bt(t),t._a=l,l.provide(hi,t),l.config.globalProperties.$pinia=t,e.forEach(o=>a.push(o)),e=[]},use(l){return this._a?a.push(l):e.push(l),this},_p:a,_a:null,_e:s,_s:new Map,state:n});return t}const di=()=>{};function Np(s,n,a,e=di){s.push(n);const t=()=>{const l=s.indexOf(n);l>-1&&(s.splice(l,1),e())};return!a&&er()&&io(t),t}function Ba(s,...n){s.slice().forEach(a=>{a(...n)})}const dm=s=>s(),Fp=Symbol(),il=Symbol();function Ol(s,n){s instanceof Map&&n instanceof Map?n.forEach((a,e)=>s.set(e,a)):s instanceof Set&&n instanceof Set&&n.forEach(s.add,s);for(const a in n){if(!n.hasOwnProperty(a))continue;const e=n[a],t=s[a];Dl(t)&&Dl(e)&&s.hasOwnProperty(a)&&!zs(e)&&!_a(e)?s[a]=Ol(t,e):s[a]=e}return s}const mm=Symbol();function jm(s){return!Dl(s)||!Object.prototype.hasOwnProperty.call(s,mm)}const{assign:da}=Object;function fm(s){return!!(zs(s)&&s.effect)}function gm(s,n,a,e){const{state:t,actions:l,getters:o}=n,p=a.state.value[s];let c;function r(){p||(a.state.value[s]=t?t():{});const i=gh(a.state.value[s]);return da(i,l,Object.keys(o||{}).reduce((u,h)=>(u[h]=_o(J(()=>{Bt(a);const d=a._s.get(s);return o[h].call(d,d)})),u),{}))}return c=mi(s,r,n,a,e,!0),c}function mi(s,n,a={},e,t,l){let o;const p=da({actions:{}},a),c={deep:!0};let r,i,u=[],h=[],d;const f=e.state.value[s];!l&&!f&&(e.state.value[s]={}),ps({});let m;function k(A){let T;r=i=!1,typeof A=="function"?(A(e.state.value[s]),T={type:ye.patchFunction,storeId:s,events:d}):(Ol(e.state.value[s],A),T={type:ye.patchObject,payload:A,storeId:s,events:d});const z=m=Symbol();vn().then(()=>{m===z&&(r=!0)}),i=!0,Ba(u,T,e.state.value[s])}const j=l?function(){const{state:T}=a,z=T?T():{};this.$patch(I=>{da(I,z)})}:di;function v(){o.stop(),u=[],h=[],e._s.delete(s)}const y=(A,T="")=>{if(Fp in A)return A[il]=T,A;const z=function(){Bt(e);const I=Array.from(arguments),D=[],W=[];function ls(es){D.push(es)}function ms(es){W.push(es)}Ba(h,{args:I,name:z[il],store:R,after:ls,onError:ms});let us;try{us=A.apply(this&&this.$id===s?this:R,I)}catch(es){throw Ba(W,es),es}return us instanceof Promise?us.then(es=>(Ba(D,es),es)).catch(es=>(Ba(W,es),Promise.reject(es))):(Ba(D,us),us)};return z[Fp]=!0,z[il]=T,z},w={_p:e,$id:s,$onAction:Np.bind(null,h),$patch:k,$reset:j,$subscribe(A,T={}){const z=Np(u,A,T.detached,()=>I()),I=o.run(()=>Ws(()=>e.state.value[s],D=>{(T.flush==="sync"?i:r)&&A({storeId:s,type:ye.direct,events:d},D)},da({},c,T)));return z},$dispose:v},R=$e(w);e._s.set(s,R);const S=(e._a&&e._a.runWithContext||dm)(()=>e._e.run(()=>(o=ro()).run(()=>n({action:y}))));for(const A in S){const T=S[A];if(zs(T)&&!fm(T)||_a(T))l||(f&&jm(T)&&(zs(T)?T.value=f[A]:Ol(T,f[A])),e.state.value[s][A]=T);else if(typeof T=="function"){const z=y(T,A);S[A]=z,p.actions[A]=T}}return da(R,S),da(ws(R),S),Object.defineProperty(R,"$state",{get:()=>e.state.value[s],set:A=>{k(T=>{da(T,A)})}}),e._p.forEach(A=>{da(R,o.run(()=>A({store:R,app:e._a,pinia:e,options:p})))}),f&&l&&a.hydrate&&a.hydrate(R.$state,f),r=!0,i=!0,R}/*! #__NO_SIDE_EFFECTS__ */function ji(s,n,a){let e;const t=typeof n=="function";e=t?a:n;function l(o,p){const c=Hr();return o=o||(c?_n(hi,null):null),o&&Bt(o),o=ui,o._s.has(s)||(t?mi(s,n,e,o):gm(s,e,o)),o._s.get(s)}return l.$id=s,l}const bm=new Set(["link","style","script","noscript"]),_m=new Set(["title","titleTemplate","script","style","noscript"]),Il=new Set(["base","meta","link","style","script","noscript"]),ym=new Set(["title","base","htmlAttrs","bodyAttrs","meta","link","style","script","noscript"]),vm=new Set(["base","title","titleTemplate","bodyAttrs","htmlAttrs","templateParams"]),wm=new Set(["key","tagPosition","tagPriority","tagDuplicateStrategy","innerHTML","textContent","processTemplateParams"]),Cm=new Set(["templateParams","htmlAttrs","bodyAttrs"]),km=new Set(["theme-color","google-site-verification","og","article","book","profile","twitter","author"]);function Nl(s,n={},a){for(const e in s){const t=s[e],l=a?`${a}:${e}`:e;typeof t=="object"&&t!==null?Nl(t,n,l):typeof t=="function"&&(n[l]=t)}return n}const fi=(()=>{if(console.createTask)return console.createTask;const s={run:n=>n()};return()=>s})();function gi(s,n,a,e){for(let t=a;t<s.length;t+=1)try{const l=e?e.run(()=>s[t](...n)):s[t](...n);if(l&&typeof l.then=="function")return Promise.resolve(l).then(()=>gi(s,n,t+1,e))}catch(l){return Promise.reject(l)}}function xm(s,n,a){if(s.length>0)return gi(s,n,0,fi(a))}function Pm(s,n,a){if(s.length>0){const e=fi(a);return Promise.all(s.map(t=>e.run(()=>t(...n))))}}function ul(s,n){for(const a of[...s])a(n)}var Sm=class{_hooks;_before;_after;_deprecatedHooks;_deprecatedMessages;constructor(){this._hooks={},this._before=void 0,this._after=void 0,this._deprecatedMessages=void 0,this._deprecatedHooks={},this.hook=this.hook.bind(this),this.callHook=this.callHook.bind(this),this.callHookWith=this.callHookWith.bind(this)}hook(s,n,a={}){if(!s||typeof n!="function")return()=>{};const e=s;let t;for(;this._deprecatedHooks[s];)t=this._deprecatedHooks[s],s=t.to;if(t&&!a.allowDeprecated){let l=t.message;l||(l=`${e} hook has been deprecated`+(t.to?`, please use ${t.to}`:"")),this._deprecatedMessages||(this._deprecatedMessages=new Set),this._deprecatedMessages.has(l)||(console.warn(l),this._deprecatedMessages.add(l))}if(!n.name)try{Object.defineProperty(n,"name",{get:()=>"_"+s.replace(/\W+/g,"_")+"_hook_cb",configurable:!0})}catch{}return this._hooks[s]=this._hooks[s]||[],this._hooks[s].push(n),()=>{n&&(this.removeHook(s,n),n=void 0)}}hookOnce(s,n){let a,e=(...t)=>(typeof a=="function"&&a(),a=void 0,e=void 0,n(...t));return a=this.hook(s,e),a}removeHook(s,n){const a=this._hooks[s];if(a){const e=a.indexOf(n);e!==-1&&a.splice(e,1),a.length===0&&(this._hooks[s]=void 0)}}clearHook(s){this._hooks[s]=void 0}deprecateHook(s,n){this._deprecatedHooks[s]=typeof n=="string"?{to:n}:n;const a=this._hooks[s]||[];this._hooks[s]=void 0;for(const e of a)this.hook(s,e)}deprecateHooks(s){for(const n in s)this.deprecateHook(n,s[n])}addHooks(s){const n=Nl(s),a=Object.keys(n).map(e=>this.hook(e,n[e]));return()=>{for(const e of a)e();a.length=0}}removeHooks(s){const n=Nl(s);for(const a in n)this.removeHook(a,n[a])}removeAllHooks(){this._hooks={}}callHook(s,...n){return this.callHookWith(xm,s,n)}callHookParallel(s,...n){return this.callHookWith(Pm,s,n)}callHookWith(s,n,a){const e=this._before||this._after?{name:n,args:a,context:{}}:void 0;this._before&&ul(this._before,e);const t=s(this._hooks[n]?[...this._hooks[n]]:[],a,n);return t instanceof Promise?t.finally(()=>{this._after&&e&&ul(this._after,e)}):(this._after&&e&&ul(this._after,e),t)}beforeEach(s){return this._before=this._before||[],this._before.push(s),()=>{if(this._before!==void 0){const n=this._before.indexOf(s);n!==-1&&this._before.splice(n,1)}}}afterEach(s){return this._after=this._after||[],this._after.push(s),()=>{if(this._after!==void 0){const n=this._after.indexOf(s);n!==-1&&this._after.splice(n,1)}}}};function Em(){return new Sm}const Tm=["name","property","http-equiv"],Am=new Set(["viewport","description","keywords","robots"]);function bi(s){const n=s.split(":");return n.length?km.has(n[1]):!1}function Fl(s){const{props:n,tag:a}=s;if(vm.has(a))return a;if(a==="link"&&n.rel==="canonical")return"canonical";if(a==="link"&&n.rel==="alternate"){if(n.hreflang)return`alternate:${n.hreflang}`;if(n.type)return`alternate:${n.type}:${n.href||""}`}if(n.charset)return"charset";if(s.tag==="meta"){for(const e of Tm)if(n[e]!==void 0){const t=n[e],l=t&&typeof t=="string"&&t.includes(":"),o=t&&Am.has(t),c=!(l||o)&&s.key?`:key:${s.key}`:"";return`${a}:${t}${c}`}}if(s.key)return`${a}:key:${s.key}`;if(n.id)return`${a}:id:${n.id}`;if(a==="link"&&n.rel==="alternate")return`alternate:${n.href||""}`;if(_m.has(a)){const e=s.textContent||s.innerHTML;if(e)return`${a}:content:${e}`}}function _i(s){const n=s._h||s._d;if(n)return n;const a=s.textContent||s.innerHTML;return a||`${s.tag}:${Object.entries(s.props).map(([e,t])=>`${e}:${String(t)}`).join(",")}`}function kt(s,n,a){typeof s==="function"&&(!a||a!=="titleTemplate"&&!(a[0]==="o"&&a[1]==="n"))&&(s=s());const t=n?n(a,s):s;if(Array.isArray(t))return t.map(l=>kt(l,n));if(t?.constructor===Object){const l={};for(const o of Object.keys(t))l[o]=kt(t[o],n,o);return l}return t}function Rm(s,n){const a=s==="style"?new Map:new Set;function e(t){if(t==null||t===void 0)return;const l=String(t).trim();if(l)if(s==="style"){const[o,...p]=l.split(":").map(c=>c?c.trim():"");o&&p.length&&a.set(o,p.join(":"))}else l.split(" ").filter(Boolean).forEach(o=>a.add(o))}return typeof n=="string"?s==="style"?n.split(";").forEach(e):e(n):Array.isArray(n)?n.forEach(t=>e(t)):n&&typeof n=="object"&&Object.entries(n).forEach(([t,l])=>{l&&l!=="false"&&(s==="style"?a.set(String(t).trim(),String(l)):e(t))}),a}function yi(s,n){if(s.props=s.props||{},!n)return s;if(s.tag==="templateParams")return s.props=n,s;const a=Il.has(s.tag)||s.tag==="htmlAttrs"||s.tag==="bodyAttrs";return Object.entries(n).forEach(([e,t])=>{if(e==="__proto__"||e==="constructor"||e==="prototype")return;if(t===null){s.props[e]=null;return}if(e==="class"||e==="style"){s.props[e]=Rm(e,t);return}if(wm.has(e)){if((e==="textContent"||e==="innerHTML")&&typeof t=="object"){let r=n.type;if(n.type||(r="application/json"),!r?.endsWith("json")&&r!=="speculationrules")return;n.type=r,s.props.type=r,s[e]=JSON.stringify(t)}else s[e]=t;return}const l=e.startsWith("data-"),o=a&&!l?e.toLowerCase():e,p=String(t),c=s.tag==="meta"&&o==="content";p==="true"||p===""?s.props[o]=l||c?p:!0:!t&&l&&p==="false"?s.props[o]="false":t!==void 0&&(s.props[o]=t)}),s}function Lm(s,n){const a=typeof n=="object"&&typeof n!="function"?n:{[s==="script"||s==="noscript"||s==="style"?"innerHTML":"textContent"]:n},e=yi({tag:s,props:{}},a);return e.key&&bm.has(e.tag)&&(e.props["data-hid"]=e._h=e.key),e.tag==="script"&&typeof e.innerHTML=="object"&&(e.innerHTML=JSON.stringify(e.innerHTML),e.props.type=e.props.type||"application/json"),Array.isArray(e.props.content)?e.props.content.map(t=>({...e,props:{...e.props,content:t}})):e}function Mm(s,n){if(!s)return[];typeof s=="function"&&(s=s());const a=(t,l)=>{for(let o=0;o<n.length;o++)l=n[o](t,l);return l};s=a(void 0,s);const e=[];return s=kt(s,a),Object.entries(s||{}).forEach(([t,l])=>{if(l!==void 0)for(const o of Array.isArray(l)?l:[l])e.push(Lm(t,o))}),e.flat()}const $p=(s,n)=>s._w===n._w?s._p-n._p:s._w-n._w,Bp={base:-10,title:10},Dm={critical:-8,high:-1,low:2},qp={meta:{"content-security-policy":-30,charset:-20,viewport:-15},link:{preconnect:20,stylesheet:60,preload:70,modulepreload:70,prefetch:90,"dns-prefetch":90,prerender:90},script:{async:30,defer:80,sync:50},style:{imported:40,sync:60}},Om=/@import/,re=s=>s===""||s===!0;function Im(s,n){if(typeof n.tagPriority=="number")return n.tagPriority;let a=100;const e=Dm[n.tagPriority]||0,t=s.resolvedOptions.disableCapoSorting?{link:{},script:{},style:{}}:qp;if(n.tag in Bp)a=Bp[n.tag];else if(n.tag==="meta"){const l=n.props["http-equiv"]==="content-security-policy"?"content-security-policy":n.props.charset?"charset":n.props.name==="viewport"?"viewport":null;l&&(a=qp.meta[l])}else if(n.tag==="link"&&n.props.rel)a=t.link[n.props.rel];else if(n.tag==="script"){const l=String(n.props.type);re(n.props.async)?a=t.script.async:n.props.src&&!re(n.props.defer)&&!re(n.props.async)&&l!=="module"&&!l.endsWith("json")||n.innerHTML&&!l.endsWith("json")?a=t.script.sync:(re(n.props.defer)&&n.props.src&&!re(n.props.async)||l==="module")&&(a=t.script.defer)}else n.tag==="style"&&(a=n.innerHTML&&Om.test(n.innerHTML)?t.style.imported:t.style.sync);return(a||100)+e}function zp(s,n){const a=typeof n=="function"?n(s):n,e=a.key||String(s.plugins.size+1);s.plugins.get(e)||(s.plugins.set(e,a),s.hooks.addHooks(a.hooks||{}))}function Nm(s={}){const n=Em();n.addHooks(s.hooks||{});const a=!s.document,e=new Map,t=new Map,l=new Set,o={_entryCount:1,plugins:t,dirty:!1,resolvedOptions:s,hooks:n,ssr:a,entries:e,headEntries(){return[...e.values()]},use:p=>zp(o,p),push(p,c){const r={...c||{}};delete r.head;const i=r._index??o._entryCount++,u={_i:i,input:p,options:r},h={_poll(d=!1){o.dirty=!0,!d&&l.add(i),n.callHook("entries:updated",o)},dispose(){e.delete(i)&&o.invalidate()},patch(d){(!r.mode||r.mode==="server"&&a||r.mode==="client"&&!a)&&(u.input=d,e.set(i,u),h._poll())}};return h.patch(p),h},async resolveTags(){const p={tagMap:new Map,tags:[],entries:[...o.entries.values()]};for(await n.callHook("entries:resolve",p);l.size;){const h=l.values().next().value;l.delete(h);const d=e.get(h);if(d){const f={tags:Mm(d.input,s.propResolvers||[]).map(m=>Object.assign(m,d.options)),entry:d};await n.callHook("entries:normalize",f),d._tags=f.tags.map((m,k)=>(m._w=Im(o,m),m._p=(d._i<<10)+k,m._d=Fl(m),m._d||(m._h=_i(m)),m))}}let c=!1;p.entries.flatMap(h=>(h._tags||[]).map(d=>({...d,props:{...d.props}}))).sort($p).reduce((h,d)=>{const f=d._d||d._h;if(!h.has(f))return h.set(f,d);const m=h.get(f);if((d?.tagDuplicateStrategy||(Cm.has(d.tag)?"merge":null)||(d.key&&d.key===m.key?"merge":null))==="merge"){const j={...m.props};Object.entries(d.props).forEach(([v,y])=>j[v]=v==="style"?new Map([...m.props.style||new Map,...y]):v==="class"?new Set([...m.props.class||new Set,...y]):y),h.set(f,{...d,props:j})}else d._p>>10===m._p>>10&&d.tag==="meta"&&bi(f)?(h.set(f,Object.assign([...Array.isArray(m)?m:[m],d],d)),c=!0):(d._w===m._w?d._p>m._p:d?._w<m?._w)&&h.set(f,d);return h},p.tagMap);const r=p.tagMap.get("title"),i=p.tagMap.get("titleTemplate");if(o._title=r?.textContent,i){const h=i?.textContent;if(o._titleTemplate=h,h){let d=typeof h=="function"?h(r?.textContent):h;typeof d=="string"&&!o.plugins.has("template-params")&&(d=d.replace("%s",r?.textContent||"")),r?d===null?p.tagMap.delete("title"):p.tagMap.set("title",{...r,textContent:d}):(i.tag="title",i.textContent=d)}}p.tags=Array.from(p.tagMap.values()),c&&(p.tags=p.tags.flat().sort($p)),await n.callHook("tags:beforeResolve",p),await n.callHook("tags:resolve",p),await n.callHook("tags:afterResolve",p);const u=[];for(const h of p.tags){const{innerHTML:d,tag:f,props:m}=h;if(ym.has(f)&&!(Object.keys(m).length===0&&!h.innerHTML&&!h.textContent)&&!(f==="meta"&&!m.content&&!m["http-equiv"]&&!m.charset)){if(f==="script"&&d){if(String(m.type).endsWith("json")){const k=typeof d=="string"?d:JSON.stringify(d);h.innerHTML=k.replace(/</g,"\\u003C")}else typeof d=="string"&&(h.innerHTML=d.replace(new RegExp(`</${f}`,"g"),`<\\/${f}`));h._d=Fl(h)}u.push(h)}}return u},invalidate(){for(const p of e.values())l.add(p._i);o.dirty=!0,n.callHook("entries:updated",o)}};return(s?.plugins||[]).forEach(p=>zp(o,p)),o.hooks.callHook("init",o),s.init?.forEach(p=>p&&o.push(p)),o}const Fm=(s,n)=>zs(n)?jh(n):n,vi="usehead";function $m(s){return{install(a){a.config.globalProperties.$unhead=s,a.config.globalProperties.$head=s,a.provide(vi,s)}}.install}function Bm(){if(Hr()){const s=_n(vi);if(s)return s}throw new Error("useHead() was called without provide context, ensure you call it through the setup() function.")}function qt(s,n={}){const a=n.head||Bm();return a.ssr?a.push(s||{},n):qm(a,s,n)}function qm(s,n,a={}){const e=ps(!1);let t;return rd(()=>{const o=e.value?{}:kt(n,Fm);t?t.patch(o):t=s.push(o,a)}),ca()&&(Dn(()=>{t.dispose()}),Ir(()=>{e.value=!0}),Or(()=>{e.value=!1})),t}const zm={created(){let s=!1;const n=ca();if(!n)return;const a=n.type;!a||!("head"in a)||(s=typeof a.head=="function"?()=>a.head.call(n.proxy):a.head,s&&qt(s))}};async function wi(s,n={}){const a=n.document||s.resolvedOptions.document;if(!a||!s.dirty)return;const e={shouldRender:!0,tags:[]};if(await s.hooks.callHook("dom:beforeRender",e),!!e.shouldRender)return s._domUpdatePromise||(s._domUpdatePromise=new Promise(async t=>{const l=new Map,o=new Promise(d=>{s.resolveTags().then(f=>{d(f.map(m=>{const k=l.get(m._d)||0,j={tag:m,id:(k?`${m._d}:${k}`:m._d)||m._h,shouldRender:!0};return m._d&&bi(m._d)&&l.set(m._d,k+1),j}))})});let p=s._dom;if(!p){p={title:a.title,elMap:new Map().set("htmlAttrs",a.documentElement).set("bodyAttrs",a.body)};for(const d of["body","head"]){const f=a[d]?.children;for(const m of f){const k=m.tagName.toLowerCase();if(!Il.has(k))continue;const j=yi({tag:k,props:{}},{innerHTML:m.innerHTML,...m.getAttributeNames().reduce((v,y)=>(v[y]=m.getAttribute(y),v),{})||{}});if(j.key=m.getAttribute("data-hid")||void 0,j._d=Fl(j)||_i(j),p.elMap.has(j._d)){let v=1,y=j._d;for(;p.elMap.has(y);)y=`${j._d}:${v++}`;p.elMap.set(y,m)}else p.elMap.set(j._d,m)}}}p.pendingSideEffects={...p.sideEffects},p.sideEffects={};function c(d,f,m){const k=`${d}:${f}`;p.sideEffects[k]=m,delete p.pendingSideEffects[k]}function r({id:d,$el:f,tag:m}){const k=m.tag.endsWith("Attrs");p.elMap.set(d,f),k||(m.textContent&&m.textContent!==f.textContent&&(f.textContent=m.textContent),m.innerHTML&&m.innerHTML!==f.innerHTML&&(f.innerHTML=m.innerHTML),c(d,"el",()=>{f?.remove(),p.elMap.delete(d)}));for(const j in m.props){if(!Object.prototype.hasOwnProperty.call(m.props,j))continue;const v=m.props[j];if(j.startsWith("on")&&typeof v=="function"){const w=f?.dataset;if(w&&w[`${j}fired`]){const R=j.slice(0,-5);v.call(f,new Event(R.substring(2)))}f.getAttribute(`data-${j}`)!==""&&((m.tag==="bodyAttrs"?a.defaultView:f).addEventListener(j.substring(2),v.bind(f)),f.setAttribute(`data-${j}`,""));continue}const y=`attr:${j}`;if(j==="class"){if(!v)continue;for(const w of v)k&&c(d,`${y}:${w}`,()=>f.classList.remove(w)),!f.classList.contains(w)&&f.classList.add(w)}else if(j==="style"){if(!v)continue;for(const[w,R]of v)c(d,`${y}:${w}`,()=>{f.style.removeProperty(w)}),f.style.setProperty(w,R)}else v!==!1&&v!==null&&(f.getAttribute(j)!==v&&f.setAttribute(j,v===!0?"":String(v)),k&&c(d,y,()=>f.removeAttribute(j)))}}const i=[],u={bodyClose:void 0,bodyOpen:void 0,head:void 0},h=await o;for(const d of h){const{tag:f,shouldRender:m,id:k}=d;if(m){if(f.tag==="title"){a.title=f.textContent,c("title","",()=>a.title=p.title);continue}d.$el=d.$el||p.elMap.get(k),d.$el?r(d):Il.has(f.tag)&&i.push(d)}}for(const d of i){const f=d.tag.tagPosition||"head";d.$el=a.createElement(d.tag.tag),r(d),u[f]=u[f]||a.createDocumentFragment(),u[f].appendChild(d.$el)}for(const d of h)await s.hooks.callHook("dom:renderTag",d,a,c);u.head&&a.head.appendChild(u.head),u.bodyOpen&&a.body.insertBefore(u.bodyOpen,a.body.firstChild),u.bodyClose&&a.body.appendChild(u.bodyClose);for(const d in p.pendingSideEffects)p.pendingSideEffects[d]();s._dom=p,await s.hooks.callHook("dom:rendered",{renders:h}),t()}).finally(()=>{s._domUpdatePromise=void 0,s.dirty=!1})),s._domUpdatePromise}function Um(s={}){const n=s.domOptions?.render||wi;s.document=s.document||(typeof window<"u"?document:void 0);const a=s.document?.head.querySelector('script[id="unhead:payload"]')?.innerHTML||!1;return Nm({...s,plugins:[...s.plugins||[],{key:"client",hooks:{"entries:updated":n}}],init:[a?JSON.parse(a):!1,...s.init||[]]})}function Hm(s,n){let a=0;return()=>{const e=++a;n(()=>{a===e&&s()})}}function Vm(s={}){const n=Um({domOptions:{render:Hm(()=>wi(n),a=>setTimeout(a,0))},...s});return n.install=$m(n),n}/*!
  * vue-router v4.5.1
  * (c) 2025 Eduardo San Martin Morote
  * @license MIT
  */const Va=typeof document<"u";function Ci(s){return typeof s=="object"||"displayName"in s||"props"in s||"__vccOpts"in s}function Gm(s){return s.__esModule||s[Symbol.toStringTag]==="Module"||s.default&&Ci(s.default)}const Ps=Object.assign;function hl(s,n){const a={};for(const e in n){const t=n[e];a[e]=Mn(t)?t.map(s):s(t)}return a}const ve=()=>{},Mn=Array.isArray,ki=/#/g,Wm=/&/g,Km=/\//g,Xm=/=/g,Ym=/\?/g,xi=/\+/g,Qm=/%5B/g,Jm=/%5D/g,Pi=/%5E/g,Zm=/%60/g,Si=/%7B/g,sj=/%7C/g,Ei=/%7D/g,nj=/%20/g;function Ro(s){return encodeURI(""+s).replace(sj,"|").replace(Qm,"[").replace(Jm,"]")}function aj(s){return Ro(s).replace(Si,"{").replace(Ei,"}").replace(Pi,"^")}function $l(s){return Ro(s).replace(xi,"%2B").replace(nj,"+").replace(ki,"%23").replace(Wm,"%26").replace(Zm,"`").replace(Si,"{").replace(Ei,"}").replace(Pi,"^")}function ej(s){return $l(s).replace(Xm,"%3D")}function tj(s){return Ro(s).replace(ki,"%23").replace(Ym,"%3F")}function lj(s){return s==null?"":tj(s).replace(Km,"%2F")}function Re(s){try{return decodeURIComponent(""+s)}catch{}return""+s}const oj=/\/$/,pj=s=>s.replace(oj,"");function dl(s,n,a="/"){let e,t={},l="",o="";const p=n.indexOf("#");let c=n.indexOf("?");return p<c&&p>=0&&(c=-1),c>-1&&(e=n.slice(0,c),l=n.slice(c+1,p>-1?p:n.length),t=s(l)),p>-1&&(e=e||n.slice(0,p),o=n.slice(p,n.length)),e=uj(e??n,a),{fullPath:e+(l&&"?")+l+o,path:e,query:t,hash:Re(o)}}function cj(s,n){const a=n.query?s(n.query):"";return n.path+(a&&"?")+a+(n.hash||"")}function Up(s,n){return!n||!s.toLowerCase().startsWith(n.toLowerCase())?s:s.slice(n.length)||"/"}function rj(s,n,a){const e=n.matched.length-1,t=a.matched.length-1;return e>-1&&e===t&&Za(n.matched[e],a.matched[t])&&Ti(n.params,a.params)&&s(n.query)===s(a.query)&&n.hash===a.hash}function Za(s,n){return(s.aliasOf||s)===(n.aliasOf||n)}function Ti(s,n){if(Object.keys(s).length!==Object.keys(n).length)return!1;for(const a in s)if(!ij(s[a],n[a]))return!1;return!0}function ij(s,n){return Mn(s)?Hp(s,n):Mn(n)?Hp(n,s):s===n}function Hp(s,n){return Mn(n)?s.length===n.length&&s.every((a,e)=>a===n[e]):s.length===1&&s[0]===n}function uj(s,n){if(s.startsWith("/"))return s;if(!s)return n;const a=n.split("/"),e=s.split("/"),t=e[e.length-1];(t===".."||t===".")&&e.push("");let l=a.length-1,o,p;for(o=0;o<e.length;o++)if(p=e[o],p!==".")if(p==="..")l>1&&l--;else break;return a.slice(0,l).join("/")+"/"+e.slice(o).join("/")}const ua={path:"/",name:void 0,params:{},query:{},hash:"",fullPath:"/",matched:[],meta:{},redirectedFrom:void 0};var Le;(function(s){s.pop="pop",s.push="push"})(Le||(Le={}));var we;(function(s){s.back="back",s.forward="forward",s.unknown=""})(we||(we={}));function hj(s){if(!s)if(Va){const n=document.querySelector("base");s=n&&n.getAttribute("href")||"/",s=s.replace(/^\w+:\/\/[^\/]+/,"")}else s="/";return s[0]!=="/"&&s[0]!=="#"&&(s="/"+s),pj(s)}const dj=/^[^#]+#/;function mj(s,n){return s.replace(dj,"#")+n}function jj(s,n){const a=document.documentElement.getBoundingClientRect(),e=s.getBoundingClientRect();return{behavior:n.behavior,left:e.left-a.left-(n.left||0),top:e.top-a.top-(n.top||0)}}const zt=()=>({left:window.scrollX,top:window.scrollY});function fj(s){let n;if("el"in s){const a=s.el,e=typeof a=="string"&&a.startsWith("#"),t=typeof a=="string"?e?document.getElementById(a.slice(1)):document.querySelector(a):a;if(!t)return;n=jj(t,s)}else n=s;"scrollBehavior"in document.documentElement.style?window.scrollTo(n):window.scrollTo(n.left!=null?n.left:window.scrollX,n.top!=null?n.top:window.scrollY)}function Vp(s,n){return(history.state?history.state.position-n:-1)+s}const Bl=new Map;function gj(s,n){Bl.set(s,n)}function bj(s){const n=Bl.get(s);return Bl.delete(s),n}let _j=()=>location.protocol+"//"+location.host;function Ai(s,n){const{pathname:a,search:e,hash:t}=n,l=s.indexOf("#");if(l>-1){let p=t.includes(s.slice(l))?s.slice(l).length:1,c=t.slice(p);return c[0]!=="/"&&(c="/"+c),Up(c,"")}return Up(a,s)+e+t}function yj(s,n,a,e){let t=[],l=[],o=null;const p=({state:h})=>{const d=Ai(s,location),f=a.value,m=n.value;let k=0;if(h){if(a.value=d,n.value=h,o&&o===f){o=null;return}k=m?h.position-m.position:0}else e(d);t.forEach(j=>{j(a.value,f,{delta:k,type:Le.pop,direction:k?k>0?we.forward:we.back:we.unknown})})};function c(){o=a.value}function r(h){t.push(h);const d=()=>{const f=t.indexOf(h);f>-1&&t.splice(f,1)};return l.push(d),d}function i(){const{history:h}=window;h.state&&h.replaceState(Ps({},h.state,{scroll:zt()}),"")}function u(){for(const h of l)h();l=[],window.removeEventListener("popstate",p),window.removeEventListener("beforeunload",i)}return window.addEventListener("popstate",p),window.addEventListener("beforeunload",i,{passive:!0}),{pauseListeners:c,listen:r,destroy:u}}function Gp(s,n,a,e=!1,t=!1){return{back:s,current:n,forward:a,replaced:e,position:window.history.length,scroll:t?zt():null}}function vj(s){const{history:n,location:a}=window,e={value:Ai(s,a)},t={value:n.state};t.value||l(e.value,{back:null,current:e.value,forward:null,position:n.length-1,replaced:!0,scroll:null},!0);function l(c,r,i){const u=s.indexOf("#"),h=u>-1?(a.host&&document.querySelector("base")?s:s.slice(u))+c:_j()+s+c;try{n[i?"replaceState":"pushState"](r,"",h),t.value=r}catch(d){console.error(d),a[i?"replace":"assign"](h)}}function o(c,r){const i=Ps({},n.state,Gp(t.value.back,c,t.value.forward,!0),r,{position:t.value.position});l(c,i,!0),e.value=c}function p(c,r){const i=Ps({},t.value,n.state,{forward:c,scroll:zt()});l(i.current,i,!0);const u=Ps({},Gp(e.value,c,null),{position:i.position+1},r);l(c,u,!1),e.value=c}return{location:e,state:t,push:p,replace:o}}function wj(s){s=hj(s);const n=vj(s),a=yj(s,n.state,n.location,n.replace);function e(l,o=!0){o||a.pauseListeners(),history.go(l)}const t=Ps({location:"",base:s,go:e,createHref:mj.bind(null,s)},n,a);return Object.defineProperty(t,"location",{enumerable:!0,get:()=>n.location.value}),Object.defineProperty(t,"state",{enumerable:!0,get:()=>n.state.value}),t}function Cj(s){return typeof s=="string"||s&&typeof s=="object"}function Ri(s){return typeof s=="string"||typeof s=="symbol"}const Li=Symbol("");var Wp;(function(s){s[s.aborted=4]="aborted",s[s.cancelled=8]="cancelled",s[s.duplicated=16]="duplicated"})(Wp||(Wp={}));function se(s,n){return Ps(new Error,{type:s,[Li]:!0},n)}function Gn(s,n){return s instanceof Error&&Li in s&&(n==null||!!(s.type&n))}const Kp="[^/]+?",kj={sensitive:!1,strict:!1,start:!0,end:!0},xj=/[.+*?^${}()[\]/\\]/g;function Pj(s,n){const a=Ps({},kj,n),e=[];let t=a.start?"^":"";const l=[];for(const r of s){const i=r.length?[]:[90];a.strict&&!r.length&&(t+="/");for(let u=0;u<r.length;u++){const h=r[u];let d=40+(a.sensitive?.25:0);if(h.type===0)u||(t+="/"),t+=h.value.replace(xj,"\\$&"),d+=40;else if(h.type===1){const{value:f,repeatable:m,optional:k,regexp:j}=h;l.push({name:f,repeatable:m,optional:k});const v=j||Kp;if(v!==Kp){d+=10;try{new RegExp(`(${v})`)}catch(w){throw new Error(`Invalid custom RegExp for param "${f}" (${v}): `+w.message)}}let y=m?`((?:${v})(?:/(?:${v}))*)`:`(${v})`;u||(y=k&&r.length<2?`(?:/${y})`:"/"+y),k&&(y+="?"),t+=y,d+=20,k&&(d+=-8),m&&(d+=-20),v===".*"&&(d+=-50)}i.push(d)}e.push(i)}if(a.strict&&a.end){const r=e.length-1;e[r][e[r].length-1]+=.7000000000000001}a.strict||(t+="/?"),a.end?t+="$":a.strict&&!t.endsWith("/")&&(t+="(?:/|$)");const o=new RegExp(t,a.sensitive?"":"i");function p(r){const i=r.match(o),u={};if(!i)return null;for(let h=1;h<i.length;h++){const d=i[h]||"",f=l[h-1];u[f.name]=d&&f.repeatable?d.split("/"):d}return u}function c(r){let i="",u=!1;for(const h of s){(!u||!i.endsWith("/"))&&(i+="/"),u=!1;for(const d of h)if(d.type===0)i+=d.value;else if(d.type===1){const{value:f,repeatable:m,optional:k}=d,j=f in r?r[f]:"";if(Mn(j)&&!m)throw new Error(`Provided param "${f}" is an array but it is not repeatable (* or + modifiers)`);const v=Mn(j)?j.join("/"):j;if(!v)if(k)h.length<2&&(i.endsWith("/")?i=i.slice(0,-1):u=!0);else throw new Error(`Missing required param "${f}"`);i+=v}}return i||"/"}return{re:o,score:e,keys:l,parse:p,stringify:c}}function Sj(s,n){let a=0;for(;a<s.length&&a<n.length;){const e=n[a]-s[a];if(e)return e;a++}return s.length<n.length?s.length===1&&s[0]===80?-1:1:s.length>n.length?n.length===1&&n[0]===80?1:-1:0}function Mi(s,n){let a=0;const e=s.score,t=n.score;for(;a<e.length&&a<t.length;){const l=Sj(e[a],t[a]);if(l)return l;a++}if(Math.abs(t.length-e.length)===1){if(Xp(e))return 1;if(Xp(t))return-1}return t.length-e.length}function Xp(s){const n=s[s.length-1];return s.length>0&&n[n.length-1]<0}const Ej={type:0,value:""},Tj=/[a-zA-Z0-9_]/;function Aj(s){if(!s)return[[]];if(s==="/")return[[Ej]];if(!s.startsWith("/"))throw new Error(`Invalid path "${s}"`);function n(d){throw new Error(`ERR (${a})/"${r}": ${d}`)}let a=0,e=a;const t=[];let l;function o(){l&&t.push(l),l=[]}let p=0,c,r="",i="";function u(){r&&(a===0?l.push({type:0,value:r}):a===1||a===2||a===3?(l.length>1&&(c==="*"||c==="+")&&n(`A repeatable param (${r}) must be alone in its segment. eg: '/:ids+.`),l.push({type:1,value:r,regexp:i,repeatable:c==="*"||c==="+",optional:c==="*"||c==="?"})):n("Invalid state to consume buffer"),r="")}function h(){r+=c}for(;p<s.length;){if(c=s[p++],c==="\\"&&a!==2){e=a,a=4;continue}switch(a){case 0:c==="/"?(r&&u(),o()):c===":"?(u(),a=1):h();break;case 4:h(),a=e;break;case 1:c==="("?a=2:Tj.test(c)?h():(u(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--);break;case 2:c===")"?i[i.length-1]=="\\"?i=i.slice(0,-1)+c:a=3:i+=c;break;case 3:u(),a=0,c!=="*"&&c!=="?"&&c!=="+"&&p--,i="";break;default:n("Unknown state");break}}return a===2&&n(`Unfinished custom RegExp for param "${r}"`),u(),o(),t}function Rj(s,n,a){const e=Pj(Aj(s.path),a),t=Ps(e,{record:s,parent:n,children:[],alias:[]});return n&&!t.record.aliasOf==!n.record.aliasOf&&n.children.push(t),t}function Lj(s,n){const a=[],e=new Map;n=Zp({strict:!1,end:!0,sensitive:!1},n);function t(u){return e.get(u)}function l(u,h,d){const f=!d,m=Qp(u);m.aliasOf=d&&d.record;const k=Zp(n,u),j=[m];if("alias"in u){const w=typeof u.alias=="string"?[u.alias]:u.alias;for(const R of w)j.push(Qp(Ps({},m,{components:d?d.record.components:m.components,path:R,aliasOf:d?d.record:m})))}let v,y;for(const w of j){const{path:R}=w;if(h&&R[0]!=="/"){const N=h.record.path,S=N[N.length-1]==="/"?"":"/";w.path=h.record.path+(R&&S+R)}if(v=Rj(w,h,k),d?d.alias.push(v):(y=y||v,y!==v&&y.alias.push(v),f&&u.name&&!Jp(v)&&o(u.name)),Di(v)&&c(v),m.children){const N=m.children;for(let S=0;S<N.length;S++)l(N[S],v,d&&d.children[S])}d=d||v}return y?()=>{o(y)}:ve}function o(u){if(Ri(u)){const h=e.get(u);h&&(e.delete(u),a.splice(a.indexOf(h),1),h.children.forEach(o),h.alias.forEach(o))}else{const h=a.indexOf(u);h>-1&&(a.splice(h,1),u.record.name&&e.delete(u.record.name),u.children.forEach(o),u.alias.forEach(o))}}function p(){return a}function c(u){const h=Oj(u,a);a.splice(h,0,u),u.record.name&&!Jp(u)&&e.set(u.record.name,u)}function r(u,h){let d,f={},m,k;if("name"in u&&u.name){if(d=e.get(u.name),!d)throw se(1,{location:u});k=d.record.name,f=Ps(Yp(h.params,d.keys.filter(y=>!y.optional).concat(d.parent?d.parent.keys.filter(y=>y.optional):[]).map(y=>y.name)),u.params&&Yp(u.params,d.keys.map(y=>y.name))),m=d.stringify(f)}else if(u.path!=null)m=u.path,d=a.find(y=>y.re.test(m)),d&&(f=d.parse(m),k=d.record.name);else{if(d=h.name?e.get(h.name):a.find(y=>y.re.test(h.path)),!d)throw se(1,{location:u,currentLocation:h});k=d.record.name,f=Ps({},h.params,u.params),m=d.stringify(f)}const j=[];let v=d;for(;v;)j.unshift(v.record),v=v.parent;return{name:k,path:m,params:f,matched:j,meta:Dj(j)}}s.forEach(u=>l(u));function i(){a.length=0,e.clear()}return{addRoute:l,resolve:r,removeRoute:o,clearRoutes:i,getRoutes:p,getRecordMatcher:t}}function Yp(s,n){const a={};for(const e of n)e in s&&(a[e]=s[e]);return a}function Qp(s){const n={path:s.path,redirect:s.redirect,name:s.name,meta:s.meta||{},aliasOf:s.aliasOf,beforeEnter:s.beforeEnter,props:Mj(s),children:s.children||[],instances:{},leaveGuards:new Set,updateGuards:new Set,enterCallbacks:{},components:"components"in s?s.components||null:s.component&&{default:s.component}};return Object.defineProperty(n,"mods",{value:{}}),n}function Mj(s){const n={},a=s.props||!1;if("component"in s)n.default=a;else for(const e in s.components)n[e]=typeof a=="object"?a[e]:a;return n}function Jp(s){for(;s;){if(s.record.aliasOf)return!0;s=s.parent}return!1}function Dj(s){return s.reduce((n,a)=>Ps(n,a.meta),{})}function Zp(s,n){const a={};for(const e in s)a[e]=e in n?n[e]:s[e];return a}function Oj(s,n){let a=0,e=n.length;for(;a!==e;){const l=a+e>>1;Mi(s,n[l])<0?e=l:a=l+1}const t=Ij(s);return t&&(e=n.lastIndexOf(t,e-1)),e}function Ij(s){let n=s;for(;n=n.parent;)if(Di(n)&&Mi(s,n)===0)return n}function Di({record:s}){return!!(s.name||s.components&&Object.keys(s.components).length||s.redirect)}function Nj(s){const n={};if(s===""||s==="?")return n;const e=(s[0]==="?"?s.slice(1):s).split("&");for(let t=0;t<e.length;++t){const l=e[t].replace(xi," "),o=l.indexOf("="),p=Re(o<0?l:l.slice(0,o)),c=o<0?null:Re(l.slice(o+1));if(p in n){let r=n[p];Mn(r)||(r=n[p]=[r]),r.push(c)}else n[p]=c}return n}function sc(s){let n="";for(let a in s){const e=s[a];if(a=ej(a),e==null){e!==void 0&&(n+=(n.length?"&":"")+a);continue}(Mn(e)?e.map(l=>l&&$l(l)):[e&&$l(e)]).forEach(l=>{l!==void 0&&(n+=(n.length?"&":"")+a,l!=null&&(n+="="+l))})}return n}function Fj(s){const n={};for(const a in s){const e=s[a];e!==void 0&&(n[a]=Mn(e)?e.map(t=>t==null?null:""+t):e==null?e:""+e)}return n}const $j=Symbol(""),nc=Symbol(""),Ut=Symbol(""),Lo=Symbol(""),ql=Symbol("");function ie(){let s=[];function n(e){return s.push(e),()=>{const t=s.indexOf(e);t>-1&&s.splice(t,1)}}function a(){s=[]}return{add:n,list:()=>s.slice(),reset:a}}function fa(s,n,a,e,t,l=o=>o()){const o=e&&(e.enterCallbacks[t]=e.enterCallbacks[t]||[]);return()=>new Promise((p,c)=>{const r=h=>{h===!1?c(se(4,{from:a,to:n})):h instanceof Error?c(h):Cj(h)?c(se(2,{from:n,to:h})):(o&&e.enterCallbacks[t]===o&&typeof h=="function"&&o.push(h),p())},i=l(()=>s.call(e&&e.instances[t],n,a,r));let u=Promise.resolve(i);s.length<3&&(u=u.then(r)),u.catch(h=>c(h))})}function ml(s,n,a,e,t=l=>l()){const l=[];for(const o of s)for(const p in o.components){let c=o.components[p];if(!(n!=="beforeRouteEnter"&&!o.instances[p]))if(Ci(c)){const i=(c.__vccOpts||c)[n];i&&l.push(fa(i,a,e,o,p,t))}else{let r=c();l.push(()=>r.then(i=>{if(!i)throw new Error(`Couldn't resolve component "${p}" at "${o.path}"`);const u=Gm(i)?i.default:i;o.mods[p]=i,o.components[p]=u;const d=(u.__vccOpts||u)[n];return d&&fa(d,a,e,o,p,t)()}))}}return l}function ac(s){const n=_n(Ut),a=_n(Lo),e=J(()=>{const c=Y(s.to);return n.resolve(c)}),t=J(()=>{const{matched:c}=e.value,{length:r}=c,i=c[r-1],u=a.matched;if(!i||!u.length)return-1;const h=u.findIndex(Za.bind(null,i));if(h>-1)return h;const d=ec(c[r-2]);return r>1&&ec(i)===d&&u[u.length-1].path!==d?u.findIndex(Za.bind(null,c[r-2])):h}),l=J(()=>t.value>-1&&Hj(a.params,e.value.params)),o=J(()=>t.value>-1&&t.value===a.matched.length-1&&Ti(a.params,e.value.params));function p(c={}){if(Uj(c)){const r=n[Y(s.replace)?"replace":"push"](Y(s.to)).catch(ve);return s.viewTransition&&typeof document<"u"&&"startViewTransition"in document&&document.startViewTransition(()=>r),r}return Promise.resolve()}return{route:e,href:J(()=>e.value.href),isActive:l,isExactActive:o,navigate:p}}function Bj(s){return s.length===1?s[0]:s}const qj=Es({name:"RouterLink",compatConfig:{MODE:3},props:{to:{type:[String,Object],required:!0},replace:Boolean,activeClass:String,exactActiveClass:String,custom:Boolean,ariaCurrentValue:{type:String,default:"page"},viewTransition:Boolean},useLink:ac,setup(s,{slots:n}){const a=$e(ac(s)),{options:e}=_n(Ut),t=J(()=>({[tc(s.activeClass,e.linkActiveClass,"router-link-active")]:a.isActive,[tc(s.exactActiveClass,e.linkExactActiveClass,"router-link-exact-active")]:a.isExactActive}));return()=>{const l=n.default&&Bj(n.default(a));return s.custom?l:Ve("a",{"aria-current":a.isExactActive?s.ariaCurrentValue:null,href:a.href,onClick:a.navigate,class:t.value},l)}}}),zj=qj;function Uj(s){if(!(s.metaKey||s.altKey||s.ctrlKey||s.shiftKey)&&!s.defaultPrevented&&!(s.button!==void 0&&s.button!==0)){if(s.currentTarget&&s.currentTarget.getAttribute){const n=s.currentTarget.getAttribute("target");if(/\b_blank\b/i.test(n))return}return s.preventDefault&&s.preventDefault(),!0}}function Hj(s,n){for(const a in n){const e=n[a],t=s[a];if(typeof e=="string"){if(e!==t)return!1}else if(!Mn(t)||t.length!==e.length||e.some((l,o)=>l!==t[o]))return!1}return!0}function ec(s){return s?s.aliasOf?s.aliasOf.path:s.path:""}const tc=(s,n,a)=>s??n??a,Vj=Es({name:"RouterView",inheritAttrs:!1,props:{name:{type:String,default:"default"},route:Object},compatConfig:{MODE:3},setup(s,{attrs:n,slots:a}){const e=_n(ql),t=J(()=>s.route||e.value),l=_n(nc,0),o=J(()=>{let r=Y(l);const{matched:i}=t.value;let u;for(;(u=i[r])&&!u.components;)r++;return r}),p=J(()=>t.value.matched[o.value]);rt(nc,J(()=>o.value+1)),rt($j,p),rt(ql,t);const c=ps();return Ws(()=>[c.value,p.value,s.name],([r,i,u],[h,d,f])=>{i&&(i.instances[u]=r,d&&d!==i&&r&&r===h&&(i.leaveGuards.size||(i.leaveGuards=d.leaveGuards),i.updateGuards.size||(i.updateGuards=d.updateGuards))),r&&i&&(!d||!Za(i,d)||!h)&&(i.enterCallbacks[u]||[]).forEach(m=>m(r))},{flush:"post"}),()=>{const r=t.value,i=s.name,u=p.value,h=u&&u.components[i];if(!h)return lc(a.default,{Component:h,route:r});const d=u.props[i],f=d?d===!0?r.params:typeof d=="function"?d(r):d:null,k=Ve(h,Ps({},f,n,{onVnodeUnmounted:j=>{j.component.isUnmounted&&(u.instances[i]=null)},ref:c}));return lc(a.default,{Component:k,route:r})||k}}});function lc(s,n){if(!s)return null;const a=s(n);return a.length===1?a[0]:a}const Gj=Vj;function Wj(s){const n=Lj(s.routes,s),a=s.parseQuery||Nj,e=s.stringifyQuery||sc,t=s.history,l=ie(),o=ie(),p=ie(),c=yo(ua);let r=ua;Va&&s.scrollBehavior&&"scrollRestoration"in history&&(history.scrollRestoration="manual");const i=hl.bind(null,F=>""+F),u=hl.bind(null,lj),h=hl.bind(null,Re);function d(F,X){let K,as;return Ri(F)?(K=n.getRecordMatcher(F),as=X):as=F,n.addRoute(as,K)}function f(F){const X=n.getRecordMatcher(F);X&&n.removeRoute(X)}function m(){return n.getRoutes().map(F=>F.record)}function k(F){return!!n.getRecordMatcher(F)}function j(F,X){if(X=Ps({},X||c.value),typeof F=="string"){const E=dl(a,F,X.path),B=n.resolve({path:E.path},X),H=t.createHref(E.fullPath);return Ps(E,B,{params:h(B.params),hash:Re(E.hash),redirectedFrom:void 0,href:H})}let K;if(F.path!=null)K=Ps({},F,{path:dl(a,F.path,X.path).path});else{const E=Ps({},F.params);for(const B in E)E[B]==null&&delete E[B];K=Ps({},F,{params:u(E)}),X.params=u(X.params)}const as=n.resolve(K,X),xs=F.hash||"";as.params=i(h(as.params));const _=cj(e,Ps({},F,{hash:aj(xs),path:as.path})),C=t.createHref(_);return Ps({fullPath:_,hash:xs,query:e===sc?Fj(F.query):F.query||{}},as,{redirectedFrom:void 0,href:C})}function v(F){return typeof F=="string"?dl(a,F,c.value.path):Ps({},F)}function y(F,X){if(r!==F)return se(8,{from:X,to:F})}function w(F){return S(F)}function R(F){return w(Ps(v(F),{replace:!0}))}function N(F){const X=F.matched[F.matched.length-1];if(X&&X.redirect){const{redirect:K}=X;let as=typeof K=="function"?K(F):K;return typeof as=="string"&&(as=as.includes("?")||as.includes("#")?as=v(as):{path:as},as.params={}),Ps({query:F.query,hash:F.hash,params:as.path!=null?{}:F.params},as)}}function S(F,X){const K=r=j(F),as=c.value,xs=F.state,_=F.force,C=F.replace===!0,E=N(K);if(E)return S(Ps(v(E),{state:typeof E=="object"?Ps({},xs,E.state):xs,force:_,replace:C}),X||K);const B=K;B.redirectedFrom=X;let H;return!_&&rj(e,as,K)&&(H=se(16,{to:B,from:as}),qs(as,as,!0,!1)),(H?Promise.resolve(H):z(B,as)).catch(U=>Gn(U)?Gn(U,2)?U:Bs(U):cs(U,B,as)).then(U=>{if(U){if(Gn(U,2))return S(Ps({replace:C},v(U.to),{state:typeof U.to=="object"?Ps({},xs,U.to.state):xs,force:_}),X||B)}else U=D(B,as,!0,C,xs);return I(B,as,U),U})}function A(F,X){const K=y(F,X);return K?Promise.reject(K):Promise.resolve()}function T(F){const X=_s.values().next().value;return X&&typeof X.runWithContext=="function"?X.runWithContext(F):F()}function z(F,X){let K;const[as,xs,_]=Kj(F,X);K=ml(as.reverse(),"beforeRouteLeave",F,X);for(const E of as)E.leaveGuards.forEach(B=>{K.push(fa(B,F,X))});const C=A.bind(null,F,X);return K.push(C),Ts(K).then(()=>{K=[];for(const E of l.list())K.push(fa(E,F,X));return K.push(C),Ts(K)}).then(()=>{K=ml(xs,"beforeRouteUpdate",F,X);for(const E of xs)E.updateGuards.forEach(B=>{K.push(fa(B,F,X))});return K.push(C),Ts(K)}).then(()=>{K=[];for(const E of _)if(E.beforeEnter)if(Mn(E.beforeEnter))for(const B of E.beforeEnter)K.push(fa(B,F,X));else K.push(fa(E.beforeEnter,F,X));return K.push(C),Ts(K)}).then(()=>(F.matched.forEach(E=>E.enterCallbacks={}),K=ml(_,"beforeRouteEnter",F,X,T),K.push(C),Ts(K))).then(()=>{K=[];for(const E of o.list())K.push(fa(E,F,X));return K.push(C),Ts(K)}).catch(E=>Gn(E,8)?E:Promise.reject(E))}function I(F,X,K){p.list().forEach(as=>T(()=>as(F,X,K)))}function D(F,X,K,as,xs){const _=y(F,X);if(_)return _;const C=X===ua,E=Va?history.state:{};K&&(as||C?t.replace(F.fullPath,Ps({scroll:C&&E&&E.scroll},xs)):t.push(F.fullPath,xs)),c.value=F,qs(F,X,K,C),Bs()}let W;function ls(){W||(W=t.listen((F,X,K)=>{if(!Ds.listening)return;const as=j(F),xs=N(as);if(xs){S(Ps(xs,{replace:!0,force:!0}),as).catch(ve);return}r=as;const _=c.value;Va&&gj(Vp(_.fullPath,K.delta),zt()),z(as,_).catch(C=>Gn(C,12)?C:Gn(C,2)?(S(Ps(v(C.to),{force:!0}),as).then(E=>{Gn(E,20)&&!K.delta&&K.type===Le.pop&&t.go(-1,!1)}).catch(ve),Promise.reject()):(K.delta&&t.go(-K.delta,!1),cs(C,as,_))).then(C=>{C=C||D(as,_,!1),C&&(K.delta&&!Gn(C,8)?t.go(-K.delta,!1):K.type===Le.pop&&Gn(C,20)&&t.go(-1,!1)),I(as,_,C)}).catch(ve)}))}let ms=ie(),us=ie(),es;function cs(F,X,K){Bs(F);const as=us.list();return as.length?as.forEach(xs=>xs(F,X,K)):console.error(F),Promise.reject(F)}function $s(){return es&&c.value!==ua?Promise.resolve():new Promise((F,X)=>{ms.add([F,X])})}function Bs(F){return es||(es=!F,ls(),ms.list().forEach(([X,K])=>F?K(F):X()),ms.reset()),F}function qs(F,X,K,as){const{scrollBehavior:xs}=s;if(!Va||!xs)return Promise.resolve();const _=!K&&bj(Vp(F.fullPath,0))||(as||!K)&&history.state&&history.state.scroll||null;return vn().then(()=>xs(F,X,_)).then(C=>C&&fj(C)).catch(C=>cs(C,F,X))}const os=F=>t.go(F);let gs;const _s=new Set,Ds={currentRoute:c,listening:!0,addRoute:d,removeRoute:f,clearRoutes:n.clearRoutes,hasRoute:k,getRoutes:m,resolve:j,options:s,push:w,replace:R,go:os,back:()=>os(-1),forward:()=>os(1),beforeEach:l.add,beforeResolve:o.add,afterEach:p.add,onError:us.add,isReady:$s,install(F){const X=this;F.component("RouterLink",zj),F.component("RouterView",Gj),F.config.globalProperties.$router=X,Object.defineProperty(F.config.globalProperties,"$route",{enumerable:!0,get:()=>Y(c)}),Va&&!gs&&c.value===ua&&(gs=!0,w(t.location).catch(xs=>{}));const K={};for(const xs in ua)Object.defineProperty(K,xs,{get:()=>c.value[xs],enumerable:!0});F.provide(Ut,X),F.provide(Lo,br(K)),F.provide(ql,c);const as=F.unmount;_s.add(F),F.unmount=function(){_s.delete(F),_s.size<1&&(r=ua,W&&W(),W=null,c.value=ua,gs=!1,es=!1),as()}}};function Ts(F){return F.reduce((X,K)=>X.then(()=>T(K)),Promise.resolve())}return Ds}function Kj(s,n){const a=[],e=[],t=[],l=Math.max(n.matched.length,s.matched.length);for(let o=0;o<l;o++){const p=n.matched[o];p&&(s.matched.find(r=>Za(r,p))?e.push(p):a.push(p));const c=s.matched[o];c&&(n.matched.find(r=>Za(r,c))||t.push(c))}return[a,e,t]}function Ge(){return _n(Ut)}function ra(s){return _n(Lo)}function Xj(s){return document.readyState==="loading"?new Promise(n=>{document.addEventListener("DOMContentLoaded",()=>n(s))}):Promise.resolve(s)}const Yj=Es({setup(s,{slots:n}){const a=ps(!1);return kn(()=>a.value=!0),()=>a.value?n.default&&n.default({}):n.placeholder&&n.placeholder({})}});function Qj(s){try{return JSON.parse(s||"{}")}catch(n){return console.error("[SSG] On state deserialization -",n,s),{}}}function Jj(s,n,a,e){const{transformState:t,registerComponents:l=!0,useHead:o=!0,rootContainer:p="#app"}={};async function c(r){const i=rm(s);let u;o&&i.use(u=Vm());const h=Wj({history:wj(n.base),...n}),{routes:d}=n;l&&i.component("ClientOnly",Yj);const f=[],j={app:i,head:u,isClient:!0,router:h,routes:d,onSSRAppRendered:()=>{},triggerOnSSRAppRendered:()=>Promise.all(f.map(R=>R())),initialState:{},transformState:t,routePath:r};await Xj(),j.initialState=t?.(window.__INITIAL_STATE__||{})||Qj(window.__INITIAL_STATE__),await a?.(j),i.use(h);let v,y=!0;h.beforeEach((R,N,S)=>{(y||v&&v===R.path)&&(y=!1,v=R.path,R.meta.state=j.initialState),S()});const w=j.initialState;return{...j,initialState:w}}return(async()=>{const{app:r,router:i}=await c();await i.isReady(),r.mount(p,!0)})(),c}/*!
  * vue-i18n v12.0.0-alpha.3
  * (c) 2016-present kazuya kawaguchi and contributors
  * Released under the MIT License.
  */const ea=typeof window<"u";let Cn,Oa;{const s=ea&&window.performance;s&&s.mark&&s.measure&&s.clearMarks&&s.clearMeasures&&(Cn=n=>{s.mark(n)},Oa=(n,a,e)=>{s.measure(n,a,e),s.clearMarks(a),s.clearMarks(e)})}const Zj=/\{([0-9a-zA-Z]+)\}/g;function Ht(s,...n){return n.length===1&&vs(n[0])&&(n=n[0]),(!n||!n.hasOwnProperty)&&(n={}),s.replace(Zj,(a,e)=>n.hasOwnProperty(e)?n[e]:"")}const Un=(s,n=!1)=>n?Symbol.for(s):Symbol(s),sf=(s,n,a)=>nf({l:s,k:n,s:a}),nf=s=>JSON.stringify(s).replace(/\u2028/g,"\\u2028").replace(/\u2029/g,"\\u2029").replace(/\u0027/g,"\\u0027"),Xs=s=>typeof s=="number"&&isFinite(s),af=s=>Mo(s)==="[object Date]",xt=s=>Mo(s)==="[object RegExp]",Vt=s=>ys(s)&&Object.keys(s).length===0,en=Object.assign,ef=Object.create,Rs=(s=null)=>ef(s);let oc;const tf=()=>oc||(oc=typeof globalThis<"u"?globalThis:typeof self<"u"?self:typeof window<"u"?window:typeof global<"u"?global:Rs());function pc(s){return s.replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&apos;")}const lf=Object.prototype.hasOwnProperty;function An(s,n){return lf.call(s,n)}const Gs=Array.isArray,Os=s=>typeof s=="function",ns=s=>typeof s=="string",Ns=s=>typeof s=="boolean",vs=s=>s!==null&&typeof s=="object",of=s=>vs(s)&&Os(s.then)&&Os(s.catch),Oi=Object.prototype.toString,Mo=s=>Oi.call(s),ys=s=>Mo(s)==="[object Object]",pf=s=>s==null?"":Gs(s)||ys(s)&&s.toString===Oi?JSON.stringify(s,null,2):String(s);function Ii(s,n=""){return s.reduce((a,e,t)=>t===0?a+e:a+n+e,"")}const cc=2;function cf(s,n=0,a=s.length){const e=s.split(/\r?\n/);let t=0;const l=[];for(let o=0;o<e.length;o++)if(t+=e[o].length+1,t>=n){for(let p=o-cc;p<=o+cc||a>t;p++){if(p<0||p>=e.length)continue;const c=p+1;l.push(`${c}${" ".repeat(3-String(c).length)}|  ${e[p]}`);const r=e[p].length;if(p===o){const i=n-(t-r)+1,u=Math.max(1,a>t?r-i:a-n);l.push("   |  "+" ".repeat(i)+"^".repeat(u))}else if(p>o){if(a>t){const i=Math.max(Math.min(a-t,r),1);l.push("   |  "+"^".repeat(i))}t+=r+1}}break}return l.join(`
`)}function ka(s,n){typeof console<"u"&&(console.warn("[intlify] "+s),n&&console.warn(n.stack))}const rc={};function rf(s){rc[s]||(rc[s]=!0,ka(s))}function Ni(){const s=new Map;return{events:s,on(a,e){const t=s.get(a);t&&t.push(e)||s.set(a,[e])},off(a,e){const t=s.get(a);t&&t.splice(t.indexOf(e)>>>0,1)},emit(a,e){(s.get(a)||[]).slice().map(t=>t(e)),(s.get("*")||[]).slice().map(t=>t(a,e))}}}const nt=s=>!vs(s)||Gs(s);function ht(s,n){if(nt(s)||nt(n))throw new Error("Invalid value");const a=[{src:s,des:n}];for(;a.length;){const{src:e,des:t}=a.pop();Object.keys(e).forEach(l=>{l!=="__proto__"&&(vs(e[l])&&!vs(t[l])&&(t[l]=Array.isArray(e[l])?[]:Rs()),nt(t[l])||nt(e[l])?t[l]=e[l]:a.push({src:e[l],des:t[l]}))})}}function uf(s,n,a){return{line:s,column:n,offset:a}}function zl(s,n,a){return{start:s,end:n}}const hs={EXPECTED_TOKEN:1,INVALID_TOKEN_IN_PLACEHOLDER:2,UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER:3,UNKNOWN_ESCAPE_SEQUENCE:4,INVALID_UNICODE_ESCAPE_SEQUENCE:5,UNBALANCED_CLOSING_BRACE:6,UNTERMINATED_CLOSING_BRACE:7,EMPTY_PLACEHOLDER:8,NOT_ALLOW_NEST_PLACEHOLDER:9,INVALID_LINKED_FORMAT:10,MUST_HAVE_MESSAGES_IN_PLURAL:11,UNEXPECTED_EMPTY_LINKED_MODIFIER:12,UNEXPECTED_EMPTY_LINKED_KEY:13,UNEXPECTED_LEXICAL_ANALYSIS:14,UNHANDLED_CODEGEN_NODE_TYPE:15,UNHANDLED_MINIFIER_NODE_TYPE:16},hf=17,df={[hs.EXPECTED_TOKEN]:"Expected token: '{0}'",[hs.INVALID_TOKEN_IN_PLACEHOLDER]:"Invalid token in placeholder: '{0}'",[hs.UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER]:"Unterminated single quote in placeholder",[hs.UNKNOWN_ESCAPE_SEQUENCE]:"Unknown escape sequence: \\{0}",[hs.INVALID_UNICODE_ESCAPE_SEQUENCE]:"Invalid unicode escape sequence: {0}",[hs.UNBALANCED_CLOSING_BRACE]:"Unbalanced closing brace",[hs.UNTERMINATED_CLOSING_BRACE]:"Unterminated closing brace",[hs.EMPTY_PLACEHOLDER]:"Empty placeholder",[hs.NOT_ALLOW_NEST_PLACEHOLDER]:"Not allowed nest placeholder",[hs.INVALID_LINKED_FORMAT]:"Invalid linked format",[hs.MUST_HAVE_MESSAGES_IN_PLURAL]:"Plural must have messages",[hs.UNEXPECTED_EMPTY_LINKED_MODIFIER]:"Unexpected empty linked modifier",[hs.UNEXPECTED_EMPTY_LINKED_KEY]:"Unexpected empty linked key",[hs.UNEXPECTED_LEXICAL_ANALYSIS]:"Unexpected lexical analysis in token: '{0}'",[hs.UNHANDLED_CODEGEN_NODE_TYPE]:"unhandled codegen node type: '{0}'",[hs.UNHANDLED_MINIFIER_NODE_TYPE]:"unhandled mimifier node type: '{0}'"};function We(s,n,a={}){const{domain:e,messages:t,args:l}=a,o=Ht((t||df)[s]||"",...l||[]),p=new SyntaxError(String(o));return p.code=s,n&&(p.location=n),p.domain=e,p}function mf(s){throw s}const jf="minifier";function Ga(s){switch(s.t=s.type,s.type){case 0:{const n=s;Ga(n.body),n.b=n.body,delete n.body;break}case 1:{const n=s,a=n.cases;for(let e=0;e<a.length;e++)Ga(a[e]);n.c=a,delete n.cases;break}case 2:{const n=s,a=n.items;for(let e=0;e<a.length;e++)Ga(a[e]);n.i=a,delete n.items,n.static&&(n.s=n.static,delete n.static);break}case 3:case 9:case 8:case 7:{const n=s;n.value&&(n.v=n.value,delete n.value);break}case 6:{const n=s;Ga(n.key),n.k=n.key,delete n.key,n.modifier&&(Ga(n.modifier),n.m=n.modifier,delete n.modifier);break}case 5:{const n=s;n.i=n.index,delete n.index;break}case 4:{const n=s;n.k=n.key,delete n.key;break}default:throw We(hs.UNHANDLED_MINIFIER_NODE_TYPE,null,{domain:jf,args:[s.type]})}delete s.type}function ff(s){const n=s.body;return n.type===2?ic(n):n.cases.forEach(a=>ic(a)),s}function ic(s){if(s.items.length===1){const n=s.items[0];(n.type===3||n.type===9)&&(s.static=n.value,delete n.value)}else{const n=[];for(let a=0;a<s.items.length;a++){const e=s.items[a];if(!(e.type===3||e.type===9)||e.value==null)break;n.push(e.value)}if(n.length===s.items.length){s.static=Ii(n);for(let a=0;a<s.items.length;a++){const e=s.items[a];(e.type===3||e.type===9)&&delete e.value}}}}const Wn=" ",gf="\r",dn=`
`,bf="\u2028",_f="\u2029";function yf(s){const n=s;let a=0,e=1,t=1,l=0;const o=S=>n[S]===gf&&n[S+1]===dn,p=S=>n[S]===dn,c=S=>n[S]===_f,r=S=>n[S]===bf,i=S=>o(S)||p(S)||c(S)||r(S),u=()=>a,h=()=>e,d=()=>t,f=()=>l,m=S=>o(S)||c(S)||r(S)?dn:n[S],k=()=>m(a),j=()=>m(a+l);function v(){return l=0,i(a)&&(e++,t=0),o(a)&&a++,a++,t++,n[a]}function y(){return o(a+l)&&l++,l++,n[a+l]}function w(){a=0,e=1,t=1,l=0}function R(S=0){l=S}function N(){const S=a+l;for(;S!==a;)v();l=0}return{index:u,line:h,column:d,peekOffset:f,charAt:m,currentChar:k,currentPeek:j,next:v,peek:y,reset:w,resetPeek:R,skipToPeek:N}}const ha=void 0,vf=".",uc="'",wf="tokenizer";function Cf(s,n={}){const a=n.location!==!1,e=yf(s),t=()=>e.index(),l=()=>uf(e.line(),e.column(),e.index()),o=l(),p=t(),c={currentType:13,offset:p,startLoc:o,endLoc:o,lastType:13,lastOffset:p,lastStartLoc:o,lastEndLoc:o,braceNest:0,inLinked:!1,text:""},r=()=>c,{onError:i}=n;function u(g,b,P,...M){const Q=r();if(b.column+=P,b.offset+=P,i){const G=a?zl(Q.startLoc,b):null,Z=We(g,G,{domain:wf,args:M});i(Z)}}function h(g,b,P){g.endLoc=l(),g.currentType=b;const M={type:b};return a&&(M.loc=zl(g.startLoc,g.endLoc)),P!=null&&(M.value=P),M}const d=g=>h(g,13);function f(g,b){return g.currentChar()===b?(g.next(),b):(u(hs.EXPECTED_TOKEN,l(),0,b),"")}function m(g){let b="";for(;g.currentPeek()===Wn||g.currentPeek()===dn;)b+=g.currentPeek(),g.peek();return b}function k(g){const b=m(g);return g.skipToPeek(),b}function j(g){if(g===ha)return!1;const b=g.charCodeAt(0);return b>=97&&b<=122||b>=65&&b<=90||b===95}function v(g){if(g===ha)return!1;const b=g.charCodeAt(0);return b>=48&&b<=57}function y(g,b){const{currentType:P}=b;if(P!==2)return!1;m(g);const M=j(g.currentPeek());return g.resetPeek(),M}function w(g,b){const{currentType:P}=b;if(P!==2)return!1;m(g);const M=g.currentPeek()==="-"?g.peek():g.currentPeek(),Q=v(M);return g.resetPeek(),Q}function R(g,b){const{currentType:P}=b;if(P!==2)return!1;m(g);const M=g.currentPeek()===uc;return g.resetPeek(),M}function N(g,b){const{currentType:P}=b;if(P!==7)return!1;m(g);const M=g.currentPeek()===".";return g.resetPeek(),M}function S(g,b){const{currentType:P}=b;if(P!==8)return!1;m(g);const M=j(g.currentPeek());return g.resetPeek(),M}function A(g,b){const{currentType:P}=b;if(!(P===7||P===11))return!1;m(g);const M=g.currentPeek()===":";return g.resetPeek(),M}function T(g,b){const{currentType:P}=b;if(P!==9)return!1;const M=()=>{const G=g.currentPeek();return G==="{"?j(g.peek()):G==="@"||G==="|"||G===":"||G==="."||G===Wn||!G?!1:G===dn?(g.peek(),M()):I(g,!1)},Q=M();return g.resetPeek(),Q}function z(g){m(g);const b=g.currentPeek()==="|";return g.resetPeek(),b}function I(g,b=!0){const P=(Q=!1,G="")=>{const Z=g.currentPeek();return Z==="{"||Z==="@"||!Z?Q:Z==="|"?!(G===Wn||G===dn):Z===Wn?(g.peek(),P(!0,Wn)):Z===dn?(g.peek(),P(!0,dn)):!0},M=P();return b&&g.resetPeek(),M}function D(g,b){const P=g.currentChar();return P===ha?ha:b(P)?(g.next(),P):null}function W(g){const b=g.charCodeAt(0);return b>=97&&b<=122||b>=65&&b<=90||b>=48&&b<=57||b===95||b===36}function ls(g){return D(g,W)}function ms(g){const b=g.charCodeAt(0);return b>=97&&b<=122||b>=65&&b<=90||b>=48&&b<=57||b===95||b===36||b===45}function us(g){return D(g,ms)}function es(g){const b=g.charCodeAt(0);return b>=48&&b<=57}function cs(g){return D(g,es)}function $s(g){const b=g.charCodeAt(0);return b>=48&&b<=57||b>=65&&b<=70||b>=97&&b<=102}function Bs(g){return D(g,$s)}function qs(g){let b="",P="";for(;b=cs(g);)P+=b;return P}function os(g){let b="";for(;;){const P=g.currentChar();if(P==="{"||P==="}"||P==="@"||P==="|"||!P)break;if(P===Wn||P===dn)if(I(g))b+=P,g.next();else{if(z(g))break;b+=P,g.next()}else b+=P,g.next()}return b}function gs(g){k(g);let b="",P="";for(;b=us(g);)P+=b;return g.currentChar()===ha&&u(hs.UNTERMINATED_CLOSING_BRACE,l(),0),P}function _s(g){k(g);let b="";return g.currentChar()==="-"?(g.next(),b+=`-${qs(g)}`):b+=qs(g),g.currentChar()===ha&&u(hs.UNTERMINATED_CLOSING_BRACE,l(),0),b}function Ds(g){return g!==uc&&g!==dn}function Ts(g){k(g),f(g,"'");let b="",P="";for(;b=D(g,Ds);)b==="\\"?P+=F(g):P+=b;const M=g.currentChar();return M===dn||M===ha?(u(hs.UNTERMINATED_SINGLE_QUOTE_IN_PLACEHOLDER,l(),0),M===dn&&(g.next(),f(g,"'")),P):(f(g,"'"),P)}function F(g){const b=g.currentChar();switch(b){case"\\":case"'":return g.next(),`\\${b}`;case"u":return X(g,b,4);case"U":return X(g,b,6);default:return u(hs.UNKNOWN_ESCAPE_SEQUENCE,l(),0,b),""}}function X(g,b,P){f(g,b);let M="";for(let Q=0;Q<P;Q++){const G=Bs(g);if(!G){u(hs.INVALID_UNICODE_ESCAPE_SEQUENCE,l(),0,`\\${b}${M}${g.currentChar()}`);break}M+=G}return`\\${b}${M}`}function K(g){return g!=="{"&&g!=="}"&&g!==Wn&&g!==dn}function as(g){k(g);let b="",P="";for(;b=D(g,K);)P+=b;return P}function xs(g){let b="",P="";for(;b=ls(g);)P+=b;return P}function _(g){const b=P=>{const M=g.currentChar();return M==="{"||M==="@"||M==="|"||M==="("||M===")"||!M||M===Wn?P:(P+=M,g.next(),b(P))};return b("")}function C(g){k(g);const b=f(g,"|");return k(g),b}function E(g,b){let P=null;switch(g.currentChar()){case"{":return b.braceNest>=1&&u(hs.NOT_ALLOW_NEST_PLACEHOLDER,l(),0),g.next(),P=h(b,2,"{"),k(g),b.braceNest++,P;case"}":return b.braceNest>0&&b.currentType===2&&u(hs.EMPTY_PLACEHOLDER,l(),0),g.next(),P=h(b,3,"}"),b.braceNest--,b.braceNest>0&&k(g),b.inLinked&&b.braceNest===0&&(b.inLinked=!1),P;case"@":return b.braceNest>0&&u(hs.UNTERMINATED_CLOSING_BRACE,l(),0),P=B(g,b)||d(b),b.braceNest=0,P;default:{let Q=!0,G=!0,Z=!0;if(z(g))return b.braceNest>0&&u(hs.UNTERMINATED_CLOSING_BRACE,l(),0),P=h(b,1,C(g)),b.braceNest=0,b.inLinked=!1,P;if(b.braceNest>0&&(b.currentType===4||b.currentType===5||b.currentType===6))return u(hs.UNTERMINATED_CLOSING_BRACE,l(),0),b.braceNest=0,H(g,b);if(Q=y(g,b))return P=h(b,4,gs(g)),k(g),P;if(G=w(g,b))return P=h(b,5,_s(g)),k(g),P;if(Z=R(g,b))return P=h(b,6,Ts(g)),k(g),P;if(!Q&&!G&&!Z)return P=h(b,12,as(g)),u(hs.INVALID_TOKEN_IN_PLACEHOLDER,l(),0,P.value),k(g),P;break}}return P}function B(g,b){const{currentType:P}=b;let M=null;const Q=g.currentChar();switch((P===7||P===8||P===11||P===9)&&(Q===dn||Q===Wn)&&u(hs.INVALID_LINKED_FORMAT,l(),0),Q){case"@":return g.next(),M=h(b,7,"@"),b.inLinked=!0,M;case".":return k(g),g.next(),h(b,8,".");case":":return k(g),g.next(),h(b,9,":");default:return z(g)?(M=h(b,1,C(g)),b.braceNest=0,b.inLinked=!1,M):N(g,b)||A(g,b)?(k(g),B(g,b)):S(g,b)?(k(g),h(b,11,xs(g))):T(g,b)?(k(g),Q==="{"?E(g,b)||M:h(b,10,_(g))):(P===7&&u(hs.INVALID_LINKED_FORMAT,l(),0),b.braceNest=0,b.inLinked=!1,H(g,b))}}function H(g,b){let P={type:13};if(b.braceNest>0)return E(g,b)||d(b);if(b.inLinked)return B(g,b)||d(b);switch(g.currentChar()){case"{":return E(g,b)||d(b);case"}":return u(hs.UNBALANCED_CLOSING_BRACE,l(),0),g.next(),h(b,3,"}");case"@":return B(g,b)||d(b);default:{if(z(g))return P=h(b,1,C(g)),b.braceNest=0,b.inLinked=!1,P;if(I(g))return h(b,0,os(g));break}}return P}function U(){const{currentType:g,offset:b,startLoc:P,endLoc:M}=c;return c.lastType=g,c.lastOffset=b,c.lastStartLoc=P,c.lastEndLoc=M,c.offset=t(),c.startLoc=l(),e.currentChar()===ha?h(c,13):H(e,c)}return{nextToken:U,currentOffset:t,currentPosition:l,context:r}}const kf="parser",xf=/(?:\\\\|\\'|\\u([0-9a-fA-F]{4})|\\U([0-9a-fA-F]{6}))/g;function Pf(s,n,a){switch(s){case"\\\\":return"\\";case"\\'":return"'";default:{const e=parseInt(n||a,16);return e<=55295||e>=57344?String.fromCodePoint(e):"�"}}}function Sf(s={}){const n=s.location!==!1,{onError:a}=s;function e(j,v,y,w,...R){const N=j.currentPosition();if(N.offset+=w,N.column+=w,a){const S=n?zl(y,N):null,A=We(v,S,{domain:kf,args:R});a(A)}}function t(j,v,y){const w={type:j};return n&&(w.start=v,w.end=v,w.loc={start:y,end:y}),w}function l(j,v,y,w){n&&(j.end=v,j.loc&&(j.loc.end=y))}function o(j,v){const y=j.context(),w=t(3,y.offset,y.startLoc);return w.value=v,l(w,j.currentOffset(),j.currentPosition()),w}function p(j,v){const y=j.context(),{lastOffset:w,lastStartLoc:R}=y,N=t(5,w,R);return N.index=parseInt(v,10),j.nextToken(),l(N,j.currentOffset(),j.currentPosition()),N}function c(j,v){const y=j.context(),{lastOffset:w,lastStartLoc:R}=y,N=t(4,w,R);return N.key=v,j.nextToken(),l(N,j.currentOffset(),j.currentPosition()),N}function r(j,v){const y=j.context(),{lastOffset:w,lastStartLoc:R}=y,N=t(9,w,R);return N.value=v.replace(xf,Pf),j.nextToken(),l(N,j.currentOffset(),j.currentPosition()),N}function i(j){const v=j.nextToken(),y=j.context(),{lastOffset:w,lastStartLoc:R}=y,N=t(8,w,R);return v.type!==11?(e(j,hs.UNEXPECTED_EMPTY_LINKED_MODIFIER,y.lastStartLoc,0),N.value="",l(N,w,R),{nextConsumeToken:v,node:N}):(v.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,y.lastStartLoc,0,Nn(v)),N.value=v.value||"",l(N,j.currentOffset(),j.currentPosition()),{node:N})}function u(j,v){const y=j.context(),w=t(7,y.offset,y.startLoc);return w.value=v,l(w,j.currentOffset(),j.currentPosition()),w}function h(j){const v=j.context(),y=t(6,v.offset,v.startLoc);let w=j.nextToken();if(w.type===8){const R=i(j);y.modifier=R.node,w=R.nextConsumeToken||j.nextToken()}switch(w.type!==9&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(w)),w=j.nextToken(),w.type===2&&(w=j.nextToken()),w.type){case 10:w.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(w)),y.key=u(j,w.value||"");break;case 4:w.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(w)),y.key=c(j,w.value||"");break;case 5:w.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(w)),y.key=p(j,w.value||"");break;case 6:w.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(w)),y.key=r(j,w.value||"");break;default:{e(j,hs.UNEXPECTED_EMPTY_LINKED_KEY,v.lastStartLoc,0);const R=j.context(),N=t(7,R.offset,R.startLoc);return N.value="",l(N,R.offset,R.startLoc),y.key=N,l(y,R.offset,R.startLoc),{nextConsumeToken:w,node:y}}}return l(y,j.currentOffset(),j.currentPosition()),{node:y}}function d(j){const v=j.context(),y=v.currentType===1?j.currentOffset():v.offset,w=v.currentType===1?v.endLoc:v.startLoc,R=t(2,y,w);R.items=[];let N=null;do{const T=N||j.nextToken();switch(N=null,T.type){case 0:T.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(T)),R.items.push(o(j,T.value||""));break;case 5:T.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(T)),R.items.push(p(j,T.value||""));break;case 4:T.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(T)),R.items.push(c(j,T.value||""));break;case 6:T.value==null&&e(j,hs.UNEXPECTED_LEXICAL_ANALYSIS,v.lastStartLoc,0,Nn(T)),R.items.push(r(j,T.value||""));break;case 7:{const z=h(j);R.items.push(z.node),N=z.nextConsumeToken||null;break}}}while(v.currentType!==13&&v.currentType!==1);const S=v.currentType===1?v.lastOffset:j.currentOffset(),A=v.currentType===1?v.lastEndLoc:j.currentPosition();return l(R,S,A),R}function f(j,v,y,w){const R=j.context();let N=w.items.length===0;const S=t(1,v,y);S.cases=[],S.cases.push(w);do{const A=d(j);N||(N=A.items.length===0),S.cases.push(A)}while(R.currentType!==13);return N&&e(j,hs.MUST_HAVE_MESSAGES_IN_PLURAL,y,0),l(S,j.currentOffset(),j.currentPosition()),S}function m(j){const v=j.context(),{offset:y,startLoc:w}=v,R=d(j);return v.currentType===13?R:f(j,y,w,R)}function k(j){const v=Cf(j,en({},s)),y=v.context(),w=t(0,y.offset,y.startLoc);return n&&w.loc&&(w.loc.source=j),w.body=m(v),s.onCacheKey&&(w.cacheKey=s.onCacheKey(j)),y.currentType!==13&&e(v,hs.UNEXPECTED_LEXICAL_ANALYSIS,y.lastStartLoc,0,j[y.offset]||""),l(w,v.currentOffset(),v.currentPosition()),w}return{parse:k}}function Nn(s){if(s.type===13)return"EOF";const n=(s.value||"").replace(/\r?\n/gu,"\\n");return n.length>10?n.slice(0,9)+"…":n}const Ef=/<\/?[\w\s="/.':;#-\/]+>/,Tf=s=>Ef.test(s);function En(s){return vs(s)&&Do(s)===0&&(An(s,"b")||An(s,"body"))}const Fi=["b","body"];function Af(s){return xa(s,Fi)}const $i=["c","cases"];function Rf(s){return xa(s,$i,[])}const Bi=["s","static"];function Lf(s){return xa(s,Bi)}const qi=["i","items"];function Mf(s){return xa(s,qi,[])}const zi=["t","type"];function Do(s){return xa(s,zi)}const Ui=["v","value"];function at(s,n){const a=xa(s,Ui);if(a!=null)return a;throw Me(n)}const Hi=["m","modifier"];function Df(s){return xa(s,Hi)}const Vi=["k","key"];function Of(s){const n=xa(s,Vi);if(n)return n;throw Me(6)}function xa(s,n,a){for(let e=0;e<n.length;e++){const t=n[e];if(An(s,t)&&s[t]!=null)return s[t]}return a}const Gi=[...Fi,...$i,...Bi,...qi,...Vi,...Hi,...Ui,...zi];function Me(s){return new Error(`unhandled node type: ${s}`)}function jl(s){return a=>If(a,s)}function If(s,n){const a=Af(n);if(a==null)throw Me(0);if(Do(a)===1){const l=Rf(a);return s.plural(l.reduce((o,p)=>[...o,hc(s,p)],[]))}else return hc(s,a)}function hc(s,n){const a=Lf(n);if(a!=null)return s.type==="text"?a:s.normalize([a]);{const e=Mf(n).reduce((t,l)=>[...t,Ul(s,l)],[]);return s.normalize(e)}}function Ul(s,n){const a=Do(n);switch(a){case 3:return at(n,a);case 9:return at(n,a);case 4:{const e=n;if(An(e,"k")&&e.k)return s.interpolate(s.named(e.k));if(An(e,"key")&&e.key)return s.interpolate(s.named(e.key));throw Me(a)}case 5:{const e=n;if(An(e,"i")&&Xs(e.i))return s.interpolate(s.list(e.i));if(An(e,"index")&&Xs(e.index))return s.interpolate(s.list(e.index));throw Me(a)}case 6:{const e=n,t=Df(e),l=Of(e);return s.linked(Ul(s,l),t?Ul(s,t):void 0,s.type)}case 7:return at(n,a);case 8:return at(n,a);default:throw new Error(`unhandled node on format message part: ${a}`)}}const Nf="Detected HTML in '{source}' message. Recommend not using HTML messages to avoid XSS.";function Ff(s,n){n&&Tf(s)&&ka(Ht(Nf,{source:s}))}const $f=s=>s;let et=Rs();function Bf(s,n={}){let a=!1;const e=n.onError||mf;n.onError=o=>{a=!0,e(o)};const l=Sf(n).parse(s);return n.optimize&&ff(l),n.mangle&&Ga(l),{ast:l,detectError:a,code:""}}function qf(s,n){if(ns(s)){const a=Ns(n.warnHtmlMessage)?n.warnHtmlMessage:!0;Ff(s,a);const t=(n.onCacheKey||$f)(s),l=et[t];if(l)return l;const{ast:o,detectError:p}=Bf(s,{...n,location:!0,mangle:!1,optimize:!1}),c=jl(o);return p?c:et[t]=c}else{if(!En(s))return ka(`the message that is resolve with key '${n.key}' is not supported for jit compilation`),(()=>s);const a=s.cacheKey;if(a){const e=et[a];return e||(et[a]=jl(s))}else return jl(s)}}let De=null;function zf(s){De=s}function Uf(s,n,a){De&&De.emit("i18n:init",{timestamp:Date.now(),i18n:s,version:n,meta:a})}const Hf=Vf("function:translate");function Vf(s){return n=>De&&De.emit(s,n)}const on={INVALID_ARGUMENT:hf,INVALID_DATE_ARGUMENT:18,INVALID_ISO_DATE_ARGUMENT:19,NOT_SUPPORT_NON_STRING_MESSAGE:20,NOT_SUPPORT_LOCALE_PROMISE_VALUE:21,NOT_SUPPORT_LOCALE_ASYNC_FUNCTION:22,NOT_SUPPORT_LOCALE_TYPE:23},Gf=24;function Jn(s){return We(s,null,{messages:Wf})}const Wf={[on.INVALID_ARGUMENT]:"Invalid arguments",[on.INVALID_DATE_ARGUMENT]:"The date provided is an invalid Date object.Make sure your Date represents a valid date.",[on.INVALID_ISO_DATE_ARGUMENT]:"The argument provided is not a valid ISO date string",[on.NOT_SUPPORT_NON_STRING_MESSAGE]:"Not support non-string message",[on.NOT_SUPPORT_LOCALE_PROMISE_VALUE]:"cannot support promise value",[on.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION]:"cannot support async function",[on.NOT_SUPPORT_LOCALE_TYPE]:"cannot support locale type"};function Oo(s,n){return n.locale!=null?dc(n.locale):dc(s.locale)}let fl;function dc(s){if(ns(s))return s;if(Os(s)){if(s.resolvedOnce&&fl!=null)return fl;if(s.constructor.name==="Function"){const n=s();if(of(n))throw Jn(on.NOT_SUPPORT_LOCALE_PROMISE_VALUE);return fl=n}else throw Jn(on.NOT_SUPPORT_LOCALE_ASYNC_FUNCTION)}else throw Jn(on.NOT_SUPPORT_LOCALE_TYPE)}function Kf(s,n,a){return[...new Set([a,...Gs(n)?n:vs(n)?Object.keys(n):ns(n)?[n]:[a]])]}function Wi(s,n,a){const e=ns(a)?a:Pt,t=s;t.__localeChainCache||(t.__localeChainCache=new Map);let l=t.__localeChainCache.get(e);if(!l){l=[];let o=[a];for(;Gs(o);)o=mc(l,o,n);const p=Gs(n)||!ys(n)?n:n.default?n.default:null;o=ns(p)?[p]:p,Gs(o)&&mc(l,o,!1),t.__localeChainCache.set(e,l)}return l}function mc(s,n,a){let e=!0;for(let t=0;t<n.length&&Ns(e);t++){const l=n[t];ns(l)&&(e=Xf(s,n[t],a))}return e}function Xf(s,n,a){let e;const t=n.split("-");do{const l=t.join("-");e=Yf(s,l,a),t.splice(-1,1)}while(t.length&&e===!0);return e}function Yf(s,n,a){let e=!1;if(!s.includes(n)&&(e=!0,n)){e=n[n.length-1]!=="!";const t=n.replace(/!/g,"");s.push(t),(Gs(a)||ys(a))&&a[t]&&(e=a[t])}return e}const Pa=[];Pa[0]={w:[0],i:[3,0],"[":[4],o:[7]};Pa[1]={w:[1],".":[2],"[":[4],o:[7]};Pa[2]={w:[2],i:[3,0],0:[3,0]};Pa[3]={i:[3,0],0:[3,0],w:[1,1],".":[2,1],"[":[4,1],o:[7,1]};Pa[4]={"'":[5,0],'"':[6,0],"[":[4,2],"]":[1,3],o:8,l:[4,0]};Pa[5]={"'":[4,0],o:8,l:[5,0]};Pa[6]={'"':[4,0],o:8,l:[6,0]};const Qf=/^\s?(?:true|false|-?[\d.]+|'[^']*'|"[^"]*")\s?$/;function Jf(s){return Qf.test(s)}function Zf(s){const n=s.charCodeAt(0),a=s.charCodeAt(s.length-1);return n===a&&(n===34||n===39)?s.slice(1,-1):s}function sg(s){if(s==null)return"o";switch(s.charCodeAt(0)){case 91:case 93:case 46:case 34:case 39:return s;case 95:case 36:case 45:return"i";case 9:case 10:case 13:case 160:case 65279:case 8232:case 8233:return"w"}return"i"}function ng(s){const n=s.trim();return s.charAt(0)==="0"&&isNaN(parseInt(s))?!1:Jf(n)?Zf(n):"*"+n}function ag(s){const n=[];let a=-1,e=0,t=0,l,o,p,c,r,i,u;const h=[];h[0]=()=>{o===void 0?o=p:o+=p},h[1]=()=>{o!==void 0&&(n.push(o),o=void 0)},h[2]=()=>{h[0](),t++},h[3]=()=>{if(t>0)t--,e=4,h[0]();else{if(t=0,o===void 0||(o=ng(o),o===!1))return!1;h[1]()}};function d(){const f=s[a+1];if(e===5&&f==="'"||e===6&&f==='"')return a++,p="\\"+f,h[0](),!0}for(;e!==null;)if(a++,l=s[a],!(l==="\\"&&d())){if(c=sg(l),u=Pa[e],r=u[c]||u.l||8,r===8||(e=r[0],r[1]!==void 0&&(i=h[r[1]],i&&(p=l,i()===!1))))return;if(e===7)return n}}const jc=new Map;function eg(s,n){return vs(s)?s[n]:null}function tg(s,n){if(!vs(s))return null;let a=jc.get(n);if(a||(a=ag(n),a&&jc.set(n,a)),!a)return null;const e=a.length;let t=s,l=0;for(;l<e;){const o=a[l];if(Gi.includes(o)&&En(t))return null;const p=t[o];if(p===void 0||Os(t))return null;t=p,l++}return t}const bn={NOT_FOUND_KEY:1,FALLBACK_TO_TRANSLATE:2,CANNOT_FORMAT_NUMBER:3,FALLBACK_TO_NUMBER_FORMAT:4,CANNOT_FORMAT_DATE:5,FALLBACK_TO_DATE_FORMAT:6,EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER:7},lg=8,og={[bn.NOT_FOUND_KEY]:"Not found '{key}' key in '{locale}' locale messages.",[bn.FALLBACK_TO_TRANSLATE]:"Fall back to translate '{key}' key with '{target}' locale.",[bn.CANNOT_FORMAT_NUMBER]:"Cannot format a number value due to not supported Intl.NumberFormat.",[bn.FALLBACK_TO_NUMBER_FORMAT]:"Fall back to number format '{key}' key with '{target}' locale.",[bn.CANNOT_FORMAT_DATE]:"Cannot format a date value due to not supported Intl.DateTimeFormat.",[bn.FALLBACK_TO_DATE_FORMAT]:"Fall back to datetime format '{key}' key with '{target}' locale.",[bn.EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER]:"This project is using Custom Message Compiler, which is an experimental feature. It may receive breaking changes or be removed in the future."};function Ia(s,...n){return Ht(og[s],...n)}const pg="12.0.0-alpha.3",Gt=-1,Pt="en-US",St="",fc=s=>`${s.charAt(0).toLocaleUpperCase()}${s.substr(1)}`;function cg(){return{upper:(s,n)=>n==="text"&&ns(s)?s.toUpperCase():n==="vnode"&&vs(s)&&"__v_isVNode"in s?s.children.toUpperCase():s,lower:(s,n)=>n==="text"&&ns(s)?s.toLowerCase():n==="vnode"&&vs(s)&&"__v_isVNode"in s?s.children.toLowerCase():s,capitalize:(s,n)=>n==="text"&&ns(s)?fc(s):n==="vnode"&&vs(s)&&"__v_isVNode"in s?fc(s.children):s}}let Ki;function rg(s){Ki=s}let Xi;function ig(s){Xi=s}let Yi;function ug(s){Yi=s}let Qi=null;const hg=s=>{Qi=s},dg=()=>Qi;let Ji=null;const gc=s=>{Ji=s},mg=()=>Ji;let bc=0;function jg(s={}){const n=Os(s.onWarn)?s.onWarn:ka,a=ns(s.version)?s.version:pg,e=ns(s.locale)||Os(s.locale)?s.locale:Pt,t=Os(e)?Pt:e,l=Gs(s.fallbackLocale)||ys(s.fallbackLocale)||ns(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:t,o=ys(s.messages)?s.messages:gl(t),p=ys(s.datetimeFormats)?s.datetimeFormats:gl(t),c=ys(s.numberFormats)?s.numberFormats:gl(t),r=en(Rs(),s.modifiers,cg()),i=s.pluralRules||Rs(),u=Os(s.missing)?s.missing:null,h=Ns(s.missingWarn)||xt(s.missingWarn)?s.missingWarn:!0,d=Ns(s.fallbackWarn)||xt(s.fallbackWarn)?s.fallbackWarn:!0,f=!!s.fallbackFormat,m=!!s.unresolving,k=Os(s.postTranslation)?s.postTranslation:null,j=ys(s.processor)?s.processor:null,v=Ns(s.warnHtmlMessage)?s.warnHtmlMessage:!0,y=!!s.escapeParameter,w=Os(s.messageCompiler)?s.messageCompiler:Ki;Os(s.messageCompiler)&&rf(Ia(bn.EXPERIMENTAL_CUSTOM_MESSAGE_COMPILER));const R=Os(s.messageResolver)?s.messageResolver:Xi||eg,N=Os(s.localeFallbacker)?s.localeFallbacker:Yi||Kf,S=vs(s.fallbackContext)?s.fallbackContext:void 0,A=s,T=vs(A.__datetimeFormatters)?A.__datetimeFormatters:new Map,z=vs(A.__numberFormatters)?A.__numberFormatters:new Map,I=vs(A.__meta)?A.__meta:{};bc++;const D={version:a,cid:bc,locale:e,fallbackLocale:l,messages:o,modifiers:r,pluralRules:i,missing:u,missingWarn:h,fallbackWarn:d,fallbackFormat:f,unresolving:m,postTranslation:k,processor:j,warnHtmlMessage:v,escapeParameter:y,messageCompiler:w,messageResolver:R,localeFallbacker:N,fallbackContext:S,onWarn:n,__meta:I};return D.datetimeFormats=p,D.numberFormats=c,D.__datetimeFormatters=T,D.__numberFormatters=z,D.__v_emitter=A.__v_emitter!=null?A.__v_emitter:void 0,Uf(D,a,I),D}const gl=s=>({[s]:Rs()});function Wt(s,n){return s instanceof RegExp?s.test(n):s}function Zi(s,n){return s instanceof RegExp?s.test(n):s}function Io(s,n,a,e,t){const{missing:l,onWarn:o}=s;{const p=s.__v_emitter;p&&p.emit("missing",{locale:a,key:n,type:t,groupId:`${t}:${n}`})}if(l!==null){const p=l(s,a,n,t);return ns(p)?p:n}else return Zi(e,n)&&o(Ia(bn.NOT_FOUND_KEY,{key:n,locale:a})),n}function ue(s,n,a){const e=s;e.__localeChainCache=new Map,s.localeFallbacker(s,a,n)}function su(s,n){return s===n?!1:s.split("-")[0]===n.split("-")[0]}function fg(s,n){const a=n.indexOf(s);if(a===-1)return!1;for(let e=a+1;e<n.length;e++)if(su(s,n[e]))return!0;return!1}const _c=typeof Intl<"u",nu={dateTimeFormat:_c&&typeof Intl.DateTimeFormat<"u",numberFormat:_c&&typeof Intl.NumberFormat<"u"};function yc(s,...n){const{datetimeFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__datetimeFormatters:p}=s;if(!nu.dateTimeFormat)return l(Ia(bn.CANNOT_FORMAT_DATE)),St;const[c,r,i,u]=Hl(...n),h=Ns(i.missingWarn)?i.missingWarn:s.missingWarn,d=Ns(i.fallbackWarn)?i.fallbackWarn:s.fallbackWarn,f=!!i.part,m=Oo(s,i),k=o(s,t,m);if(!ns(c)||c==="")return new Intl.DateTimeFormat(m,u).format(r);let j={},v,y=null,w=m,R=null;const N="datetime format";for(let T=0;T<k.length;T++){if(v=R=k[T],m!==v&&Wt(d,c)&&l(Ia(bn.FALLBACK_TO_DATE_FORMAT,{key:c,target:v})),m!==v){const z=s.__v_emitter;z&&z.emit("fallback",{type:N,key:c,from:w,to:R,groupId:`${N}:${c}`})}if(j=a[v]||{},y=j[c],ys(y))break;Io(s,c,v,h,N),w=R}if(!ys(y)||!ns(v))return e?Gt:c;let S=`${v}__${c}`;Vt(u)||(S=`${S}__${JSON.stringify(u)}`);let A=p.get(S);return A||(A=new Intl.DateTimeFormat(v,en({},y,u)),p.set(S,A)),f?A.formatToParts(r):A.format(r)}const au=["localeMatcher","weekday","era","year","month","day","hour","minute","second","timeZoneName","formatMatcher","hour12","timeZone","dateStyle","timeStyle","calendar","dayPeriod","numberingSystem","hourCycle","fractionalSecondDigits"];function Hl(...s){const[n,a,e,t]=s,l=Rs();let o=Rs(),p;if(ns(n)){const c=n.match(/(\d{4}-\d{2}-\d{2})(T|\s)?(.*)/);if(!c)throw Jn(on.INVALID_ISO_DATE_ARGUMENT);const r=c[3]?c[3].trim().startsWith("T")?`${c[1].trim()}${c[3].trim()}`:`${c[1].trim()}T${c[3].trim()}`:c[1].trim();p=new Date(r);try{p.toISOString()}catch{throw Jn(on.INVALID_ISO_DATE_ARGUMENT)}}else if(af(n)){if(isNaN(n.getTime()))throw Jn(on.INVALID_DATE_ARGUMENT);p=n}else if(Xs(n))p=n;else throw Jn(on.INVALID_ARGUMENT);return ns(a)?l.key=a:ys(a)&&Object.keys(a).forEach(c=>{au.includes(c)?o[c]=a[c]:l[c]=a[c]}),ns(e)?l.locale=e:ys(e)&&(o=e),ys(t)&&(o=t),[l.key||"",p,l,o]}function vc(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__datetimeFormatters.has(l)&&e.__datetimeFormatters.delete(l)}}function wc(s,...n){const{numberFormats:a,unresolving:e,fallbackLocale:t,onWarn:l,localeFallbacker:o}=s,{__numberFormatters:p}=s;if(!nu.numberFormat)return l(Ia(bn.CANNOT_FORMAT_NUMBER)),St;const[c,r,i,u]=Vl(...n),h=Ns(i.missingWarn)?i.missingWarn:s.missingWarn,d=Ns(i.fallbackWarn)?i.fallbackWarn:s.fallbackWarn,f=!!i.part,m=Oo(s,i),k=o(s,t,m);if(!ns(c)||c==="")return new Intl.NumberFormat(m,u).format(r);let j={},v,y=null,w=m,R=null;const N="number format";for(let T=0;T<k.length;T++){if(v=R=k[T],m!==v&&Wt(d,c)&&l(Ia(bn.FALLBACK_TO_NUMBER_FORMAT,{key:c,target:v})),m!==v){const z=s.__v_emitter;z&&z.emit("fallback",{type:N,key:c,from:w,to:R,groupId:`${N}:${c}`})}if(j=a[v]||{},y=j[c],ys(y))break;Io(s,c,v,h,N),w=R}if(!ys(y)||!ns(v))return e?Gt:c;let S=`${v}__${c}`;Vt(u)||(S=`${S}__${JSON.stringify(u)}`);let A=p.get(S);return A||(A=new Intl.NumberFormat(v,en({},y,u)),p.set(S,A)),f?A.formatToParts(r):A.format(r)}const eu=["localeMatcher","style","currency","currencyDisplay","currencySign","useGrouping","minimumIntegerDigits","minimumFractionDigits","maximumFractionDigits","minimumSignificantDigits","maximumSignificantDigits","compactDisplay","notation","signDisplay","unit","unitDisplay","roundingMode","roundingPriority","roundingIncrement","trailingZeroDisplay"];function Vl(...s){const[n,a,e,t]=s,l=Rs();let o=Rs();if(!Xs(n))throw Jn(on.INVALID_ARGUMENT);const p=n;return ns(a)?l.key=a:ys(a)&&Object.keys(a).forEach(c=>{eu.includes(c)?o[c]=a[c]:l[c]=a[c]}),ns(e)?l.locale=e:ys(e)&&(o=e),ys(t)&&(o=t),[l.key||"",p,l,o]}function Cc(s,n,a){const e=s;for(const t in a){const l=`${n}__${t}`;e.__numberFormatters.has(l)&&e.__numberFormatters.delete(l)}}const gg=s=>s,bg=s=>"",_g="text",yg=s=>s.length===0?"":Ii(s),vg=pf;function kc(s,n){return s=Math.abs(s),n===2?s?s>1?1:0:1:s?Math.min(s,2):0}function wg(s){const n=Xs(s.pluralIndex)?s.pluralIndex:-1;return s.named&&(Xs(s.named.count)||Xs(s.named.n))?Xs(s.named.count)?s.named.count:Xs(s.named.n)?s.named.n:n:n}function Cg(s,n){n.count||(n.count=s),n.n||(n.n=s)}function kg(s={}){const n=s.locale,a=wg(s),e=vs(s.pluralRules)&&ns(n)&&Os(s.pluralRules[n])?s.pluralRules[n]:kc,t=vs(s.pluralRules)&&ns(n)&&Os(s.pluralRules[n])?kc:void 0,l=j=>j[e(a,j.length,t)],o=s.list||[],p=j=>o[j],c=s.named||Rs();Xs(s.pluralIndex)&&Cg(a,c);const r=j=>c[j];function i(j,v){const y=Os(s.messages)?s.messages(j,!!v):vs(s.messages)?s.messages[j]:!1;return y||(s.parent?s.parent.message(j):bg)}const u=j=>s.modifiers?s.modifiers[j]:gg,h=ys(s.processor)&&Os(s.processor.normalize)?s.processor.normalize:yg,d=ys(s.processor)&&Os(s.processor.interpolate)?s.processor.interpolate:vg,f=ys(s.processor)&&ns(s.processor.type)?s.processor.type:_g,k={list:p,named:r,plural:l,linked:(j,...v)=>{const[y,w]=v;let R="text",N="";v.length===1?vs(y)?(N=y.modifier||N,R=y.type||R):ns(y)&&(N=y||N):v.length===2&&(ns(y)&&(N=y||N),ns(w)&&(R=w||R));const S=i(j,!0)(k),A=R==="vnode"&&Gs(S)&&N?S[0]:S;return N?u(N)(A,R):A},message:i,type:f,interpolate:d,normalize:h,values:en(Rs(),o,c)};return k}const xc=()=>"",Pn=s=>Os(s);function Pc(s,...n){const{fallbackFormat:a,postTranslation:e,unresolving:t,messageCompiler:l,fallbackLocale:o,messages:p}=s,[c,r]=Gl(...n),i=Ns(r.missingWarn)?r.missingWarn:s.missingWarn,u=Ns(r.fallbackWarn)?r.fallbackWarn:s.fallbackWarn,h=Ns(r.escapeParameter)?r.escapeParameter:s.escapeParameter,d=!!r.resolvedMessage,f=ns(r.default)||Ns(r.default)?Ns(r.default)?l?c:()=>c:r.default:a?l?c:()=>c:null,m=a||f!=null&&(ns(f)||Os(f)),k=Oo(s,r);h&&xg(r);let[j,v,y]=d?[c,k,p[k]||Rs()]:tu(s,c,k,o,u,i),w=j,R=c;if(!d&&!(ns(w)||En(w)||Pn(w))&&m&&(w=f,R=w),!d&&(!(ns(w)||En(w)||Pn(w))||!ns(v)))return t?Gt:c;if(ns(w)&&s.messageCompiler==null)return ka(`The message format compilation is not supported in this build. Because message compiler isn't included. You need to pre-compilation all message format. So translate function return '${c}'.`),c;let N=!1;const S=()=>{N=!0},A=Pn(w)?w:lu(s,c,v,w,R,S);if(N)return w;const T=Tg(s,v,y,r),z=kg(T),I=Pg(s,A,z),D=e?e(I,c):I;{const W={timestamp:Date.now(),key:ns(c)?c:Pn(w)?w.key:"",locale:v||(Pn(w)?w.locale:""),format:ns(w)?w:Pn(w)?w.source:"",message:D};W.meta=en({},s.__meta,dg()||{}),Hf(W)}return D}function xg(s){Gs(s.list)?s.list=s.list.map(n=>ns(n)?pc(n):n):vs(s.named)&&Object.keys(s.named).forEach(n=>{ns(s.named[n])&&(s.named[n]=pc(s.named[n]))})}function tu(s,n,a,e,t,l){const{messages:o,onWarn:p,messageResolver:c,localeFallbacker:r}=s,i=r(s,e,a);let u=Rs(),h,d=null,f=a,m=null;const k="translate";for(let j=0;j<i.length;j++){if(h=m=i[j],a!==h&&!su(a,h)&&Wt(t,n)&&p(Ia(bn.FALLBACK_TO_TRANSLATE,{key:n,target:h})),a!==h){const R=s.__v_emitter;R&&R.emit("fallback",{type:k,key:n,from:f,to:m,groupId:`${k}:${n}`})}u=o[h]||Rs();let v=null,y,w;if(ea&&(v=window.performance.now(),y="intlify-message-resolve-start",w="intlify-message-resolve-end",Cn&&Cn(y)),(d=c(u,n))===null&&(d=u[n]),ea){const R=window.performance.now(),N=s.__v_emitter;N&&v&&d&&N.emit("message-resolve",{type:"message-resolve",key:n,message:d,time:R-v,groupId:`${k}:${n}`}),y&&w&&Cn&&Oa&&(Cn(w),Oa("intlify message resolve",y,w))}if(ns(d)||En(d)||Pn(d))break;if(!fg(h,i)){const R=Io(s,n,h,l,k);R!==n&&(d=R)}f=m}return[d,h,u]}function lu(s,n,a,e,t,l){const{messageCompiler:o,warnHtmlMessage:p}=s;if(Pn(e)){const h=e;return h.locale=h.locale||a,h.key=h.key||n,h}if(o==null){const h=(()=>e);return h.locale=a,h.key=n,h}let c=null,r,i;ea&&(c=window.performance.now(),r="intlify-message-compilation-start",i="intlify-message-compilation-end",Cn&&Cn(r));const u=o(e,Sg(s,a,t,e,p,l));if(ea){const h=window.performance.now(),d=s.__v_emitter;d&&c&&d.emit("message-compilation",{type:"message-compilation",message:e,time:h-c,groupId:`translate:${n}`}),r&&i&&Cn&&Oa&&(Cn(i),Oa("intlify message compilation",r,i))}return u.locale=a,u.key=n,u.source=e,u}function Pg(s,n,a){let e=null,t,l;ea&&(e=window.performance.now(),t="intlify-message-evaluation-start",l="intlify-message-evaluation-end",Cn&&Cn(t));const o=n(a);if(ea){const p=window.performance.now(),c=s.__v_emitter;c&&e&&c.emit("message-evaluation",{type:"message-evaluation",value:o,time:p-e,groupId:`translate:${n.key}`}),t&&l&&Cn&&Oa&&(Cn(l),Oa("intlify message evaluation",t,l))}return o}function Gl(...s){const[n,a,e]=s,t=Rs();if(!ns(n)&&!Xs(n)&&!Pn(n)&&!En(n))throw Jn(on.INVALID_ARGUMENT);const l=Xs(n)?String(n):(Pn(n),n);return Xs(a)?t.plural=a:ns(a)?t.default=a:ys(a)&&!Vt(a)?t.named=a:Gs(a)&&(t.list=a),Xs(e)?t.plural=e:ns(e)?t.default=e:ys(e)&&en(t,e),[l,t]}function Sg(s,n,a,e,t,l){return{locale:n,key:a,warnHtmlMessage:t,onError:o=>{l&&l(o);{const p=Eg(e),c=`Message compilation error: ${o.message}`,r=o.location&&p&&cf(p,o.location.start.offset,o.location.end.offset),i=s.__v_emitter;i&&p&&i.emit("compile-error",{message:p,error:o.message,start:o.location&&o.location.start.offset,end:o.location&&o.location.end.offset,groupId:`translate:${a}`}),console.error(r?`${c}
${r}`:c)}},onCacheKey:o=>sf(n,a,o)}}function Eg(s){if(ns(s))return s;if(s.loc&&s.loc.source)return s.loc.source}function Tg(s,n,a,e){const{modifiers:t,pluralRules:l,messageResolver:o,fallbackLocale:p,fallbackWarn:c,missingWarn:r,fallbackContext:i}=s,h={locale:n,modifiers:t,pluralRules:l,messages:(d,f)=>{let m=o(a,d);if(m==null&&(i||f)){const[,,k]=tu(i||s,d,n,p,c,r);m=o(k,d)}if(ns(m)||En(m)){let k=!1;const v=lu(s,d,n,m,d,()=>{k=!0});return k?xc:v}else return Pn(m)?m:xc}};return s.processor&&(h.processor=s.processor),e.list&&(h.list=e.list),e.named&&(h.named=e.named),Xs(e.plural)&&(h.pluralIndex=e.plural),h}const Ag="12.0.0-alpha.3";function Rg(){console.info(`You are running a development build of vue-i18n.
Make sure to use the production build (*.prod.js) when deploying for production.`)}const Ks={UNEXPECTED_RETURN_TYPE:Gf,INVALID_ARGUMENT:25,MUST_BE_CALL_SETUP_TOP:26,NOT_INSTALLED:27,REQUIRED_VALUE:28,INVALID_VALUE:29,CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN:30,NOT_INSTALLED_WITH_PROVIDE:31,UNEXPECTED_ERROR:32,DUPLICATE_USE_I18N_CALLING:33};function zn(s,...n){return We(s,null,{messages:Lg,args:n})}const Lg={[Ks.UNEXPECTED_RETURN_TYPE]:"Unexpected return type in composer",[Ks.INVALID_ARGUMENT]:"Invalid argument",[Ks.MUST_BE_CALL_SETUP_TOP]:"Must be called at the top of a `setup` function",[Ks.NOT_INSTALLED]:"Need to install with `app.use` function",[Ks.UNEXPECTED_ERROR]:"Unexpected error",[Ks.REQUIRED_VALUE]:"Required in value: {0}",[Ks.INVALID_VALUE]:"Invalid value",[Ks.CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN]:"Cannot setup vue-devtools plugin",[Ks.NOT_INSTALLED_WITH_PROVIDE]:"Need to install with `provide` function",[Ks.DUPLICATE_USE_I18N_CALLING]:"Duplicate local-scope `useI18n` call detected. Call `useI18n` only once per component."},Wl=Un("__translateVNode"),Kl=Un("__datetimeParts"),Xl=Un("__numberParts"),Oe=Un("__enableEmitter"),Yl=Un("__disableEmitter"),Mg=Un("__setPluralRules"),Dg=Un("__injectWithOption"),Ql=Un("__dispose"),Qa={FALLBACK_TO_ROOT:lg,NOT_FOUND_PARENT_SCOPE:9,IGNORE_OBJ_FLATTEN:10},Og={[Qa.FALLBACK_TO_ROOT]:"Fall back to {type} '{key}' with root locale.",[Qa.NOT_FOUND_PARENT_SCOPE]:"Not found parent scope. use the global scope.",[Qa.IGNORE_OBJ_FLATTEN]:"Ignore object flatten: '{key}' key has an string value"};function No(s,...n){return Ht(Og[s],...n)}function Ie(s){if(!vs(s)||En(s))return s;for(const n in s)if(An(s,n))if(!n.includes("."))vs(s[n])&&Ie(s[n]);else{const a=n.split("."),e=a.length-1;let t=s,l=!1;for(let o=0;o<e;o++){if(a[o]==="__proto__")throw new Error(`unsafe key: ${a[o]}`);if(a[o]in t||(t[a[o]]=Rs()),!vs(t[a[o]])){ka(No(Qa.IGNORE_OBJ_FLATTEN,{key:a[o]})),l=!0;break}t=t[a[o]]}if(l||(En(t)?Gi.includes(a[e])||delete s[n]:(t[a[e]]=s[n],delete s[n])),!En(t)){const o=t[a[e]];vs(o)&&Ie(o)}}return s}function ou(s,n){const{messages:a,__i18n:e,messageResolver:t,flatJson:l}=n,o=ys(a)?a:Gs(e)?Rs():{[s]:Rs()};if(Gs(e)&&e.forEach(p=>{if("locale"in p&&"resource"in p){const{locale:c,resource:r}=p;c?(o[c]=o[c]||Rs(),ht(r,o[c])):ht(r,o)}else ns(p)&&ht(JSON.parse(p),o)}),t==null&&l)for(const p in o)An(o,p)&&Ie(o[p]);return o}function pu(s){return s.type}function Ig(s,n,a){let e=vs(n.messages)?n.messages:Rs();"__i18nGlobal"in a&&(e=ou(s.locale.value,{messages:e,__i18n:a.__i18nGlobal}));const t=Object.keys(e);t.length&&t.forEach(l=>{s.mergeLocaleMessage(l,e[l])});{if(vs(n.datetimeFormats)){const l=Object.keys(n.datetimeFormats);l.length&&l.forEach(o=>{s.mergeDateTimeFormat(o,n.datetimeFormats[o])})}if(vs(n.numberFormats)){const l=Object.keys(n.numberFormats);l.length&&l.forEach(o=>{s.mergeNumberFormat(o,n.numberFormats[o])})}}}function Sc(s){return fs(Ue,null,s,0)}const Ec="__INTLIFY_META__",Tc=()=>[],Ng=()=>!1;let Ac=0;function Rc(s){return((n,a,e,t)=>s(a,e,ca()||void 0,t))}const Fg=()=>{const s=ca();let n=null;return s&&(n=pu(s)[Ec])?{[Ec]:n}:null};function cu(s={}){const{__root:n,__injectWithOption:a}=s,e=n===void 0,t=s.flatJson,l=ea?ps:yo;let o=Ns(s.inheritLocale)?s.inheritLocale:!0;const p=l(n&&o?n.locale.value:ns(s.locale)?s.locale:Pt),c=l(n&&o?n.fallbackLocale.value:ns(s.fallbackLocale)||Gs(s.fallbackLocale)||ys(s.fallbackLocale)||s.fallbackLocale===!1?s.fallbackLocale:p.value),r=l(ou(p.value,s)),i=l(ys(s.datetimeFormats)?s.datetimeFormats:{[p.value]:{}}),u=l(ys(s.numberFormats)?s.numberFormats:{[p.value]:{}});let h=n?n.missingWarn:Ns(s.missingWarn)||xt(s.missingWarn)?s.missingWarn:!0,d=n?n.fallbackWarn:Ns(s.fallbackWarn)||xt(s.fallbackWarn)?s.fallbackWarn:!0,f=n?n.fallbackRoot:Ns(s.fallbackRoot)?s.fallbackRoot:!0,m=!!s.fallbackFormat,k=Os(s.missing)?s.missing:null,j=Os(s.missing)?Rc(s.missing):null,v=Os(s.postTranslation)?s.postTranslation:null,y=n?n.warnHtmlMessage:Ns(s.warnHtmlMessage)?s.warnHtmlMessage:!0,w=!!s.escapeParameter;const R=n?n.modifiers:ys(s.modifiers)?s.modifiers:{};let N=s.pluralRules||n&&n.pluralRules,S;S=(()=>{e&&gc(null);const L={version:Ag,locale:p.value,fallbackLocale:c.value,messages:r.value,modifiers:R,pluralRules:N,missing:j===null?void 0:j,missingWarn:h,fallbackWarn:d,fallbackFormat:m,unresolving:!0,postTranslation:v===null?void 0:v,warnHtmlMessage:y,escapeParameter:w,messageResolver:s.messageResolver,messageCompiler:s.messageCompiler,__meta:{framework:"vue"}};L.datetimeFormats=i.value,L.numberFormats=u.value,L.__datetimeFormatters=ys(S)?S.__datetimeFormatters:void 0,L.__numberFormatters=ys(S)?S.__numberFormatters:void 0,L.__v_emitter=ys(S)?S.__v_emitter:void 0;const q=jg(L);return e&&gc(q),q})(),ue(S,p.value,c.value);function T(){return[p.value,c.value,r.value,i.value,u.value]}const z=J({get:()=>p.value,set:L=>{S.locale=L,p.value=L}}),I=J({get:()=>c.value,set:L=>{S.fallbackLocale=L,c.value=L,ue(S,p.value,L)}}),D=J(()=>r.value),W=J(()=>Object.keys(r.value).sort()),ls=J(()=>i.value),ms=J(()=>u.value);function us(){return Os(v)?v:null}function es(L){v=L,S.postTranslation=L}function cs(){return k}function $s(L){L!==null&&(j=Rc(L)),k=L,S.missing=j}function Bs(L,q){return L!=="translate"||!q.resolvedMessage}const qs=(L,q,ts,bs,Vs,un)=>{T();let Js;try{e||(S.fallbackContext=n?mg():void 0),Js=L(S)}finally{e||(S.fallbackContext=void 0)}if(ts!=="translate exists"&&Xs(Js)&&Js===Gt||ts==="translate exists"&&!Js){const[yn,Ke]=q();if(n&&ns(yn)&&Bs(ts,Ke)){f&&(Wt(d,yn)||Zi(h,yn))&&ka(No(Qa.FALLBACK_TO_ROOT,{key:yn,type:ts}));{const{__v_emitter:tn}=S;tn&&f&&tn.emit("fallback",{type:ts,key:yn,to:"global",groupId:`${ts}:${yn}`})}}return n&&f?bs(n):Vs(yn)}else{if(un(Js))return Js;throw zn(Ks.UNEXPECTED_RETURN_TYPE)}};function os(...L){return qs(q=>Reflect.apply(Pc,null,[q,...L]),()=>Gl(...L),"translate",q=>Reflect.apply(q.t,q,[...L]),q=>q,q=>ns(q))}function gs(...L){const[q,ts,bs]=L;if(bs&&!vs(bs))throw zn(Ks.INVALID_ARGUMENT);return os(q,ts,en({resolvedMessage:!0},bs||{}))}function _s(...L){return qs(q=>Reflect.apply(yc,null,[q,...L]),()=>Hl(...L),"datetime format",q=>Reflect.apply(q.d,q,[...L]),()=>St,q=>ns(q)||Gs(q))}function Ds(...L){return qs(q=>Reflect.apply(wc,null,[q,...L]),()=>Vl(...L),"number format",q=>Reflect.apply(q.n,q,[...L]),()=>St,q=>ns(q)||Gs(q))}function Ts(L){return L.map(q=>ns(q)||Xs(q)||Ns(q)?Sc(String(q)):q)}const X={normalize:Ts,interpolate:L=>L,type:"vnode"};function K(...L){return qs(q=>{let ts;const bs=q;try{bs.processor=X,ts=Reflect.apply(Pc,null,[bs,...L])}finally{bs.processor=null}return ts},()=>Gl(...L),"translate",q=>q[Wl](...L),q=>[Sc(q)],q=>Gs(q))}function as(...L){return qs(q=>Reflect.apply(wc,null,[q,...L]),()=>Vl(...L),"number format",q=>q[Xl](...L),Tc,q=>ns(q)||Gs(q))}function xs(...L){return qs(q=>Reflect.apply(yc,null,[q,...L]),()=>Hl(...L),"datetime format",q=>q[Kl](...L),Tc,q=>ns(q)||Gs(q))}function _(L){N=L,S.pluralRules=N}function C(L,q){return qs(()=>{if(!L)return!1;const ts=ns(q)?q:p.value,bs=H(ts),Vs=S.messageResolver(bs,L);return En(Vs)||Pn(Vs)||ns(Vs)},()=>[L],"translate exists",ts=>Reflect.apply(ts.te,ts,[L,q]),Ng,ts=>Ns(ts))}function E(L){let q=null;const ts=Wi(S,c.value,p.value);for(let bs=0;bs<ts.length;bs++){const Vs=r.value[ts[bs]]||{},un=S.messageResolver(Vs,L);if(un!=null){q=un;break}}return q}function B(L){const q=E(L);return q??(n?n.tm(L)||{}:{})}function H(L){return r.value[L]||{}}function U(L,q){if(t){const ts={[L]:q};for(const bs in ts)An(ts,bs)&&Ie(ts[bs]);q=ts[L]}r.value[L]=q,S.messages=r.value}function g(L,q){r.value[L]=r.value[L]||{};const ts={[L]:q};if(t)for(const bs in ts)An(ts,bs)&&Ie(ts[bs]);q=ts[L],ht(q,r.value[L]),S.messages=r.value}function b(L){return i.value[L]||{}}function P(L,q){i.value[L]=q,S.datetimeFormats=i.value,vc(S,L,q)}function M(L,q){i.value[L]=en(i.value[L]||{},q),S.datetimeFormats=i.value,vc(S,L,q)}function Q(L){return u.value[L]||{}}function G(L,q){u.value[L]=q,S.numberFormats=u.value,Cc(S,L,q)}function Z(L,q){u.value[L]=en(u.value[L]||{},q),S.numberFormats=u.value,Cc(S,L,q)}Ac++,n&&ea&&(Ws(n.locale,L=>{o&&(p.value=L,S.locale=L,ue(S,p.value,c.value))}),Ws(n.fallbackLocale,L=>{o&&(c.value=L,S.fallbackLocale=L,ue(S,p.value,c.value))}));const ss={id:Ac,locale:z,fallbackLocale:I,get inheritLocale(){return o},set inheritLocale(L){o=L,L&&n&&(p.value=n.locale.value,c.value=n.fallbackLocale.value,ue(S,p.value,c.value))},availableLocales:W,messages:D,get modifiers(){return R},get pluralRules(){return N||{}},get isGlobal(){return e},get missingWarn(){return h},set missingWarn(L){h=L,S.missingWarn=h},get fallbackWarn(){return d},set fallbackWarn(L){d=L,S.fallbackWarn=d},get fallbackRoot(){return f},set fallbackRoot(L){f=L},get fallbackFormat(){return m},set fallbackFormat(L){m=L,S.fallbackFormat=m},get warnHtmlMessage(){return y},set warnHtmlMessage(L){y=L,S.warnHtmlMessage=L},get escapeParameter(){return w},set escapeParameter(L){w=L,S.escapeParameter=L},t:os,getLocaleMessage:H,setLocaleMessage:U,mergeLocaleMessage:g,getPostTranslationHandler:us,setPostTranslationHandler:es,getMissingHandler:cs,setMissingHandler:$s,[Mg]:_};return ss.datetimeFormats=ls,ss.numberFormats=ms,ss.rt=gs,ss.te=C,ss.tm=B,ss.d=_s,ss.n=Ds,ss.getDateTimeFormat=b,ss.setDateTimeFormat=P,ss.mergeDateTimeFormat=M,ss.getNumberFormat=Q,ss.setNumberFormat=G,ss.mergeNumberFormat=Z,ss[Dg]=a,ss[Wl]=K,ss[Kl]=xs,ss[Xl]=as,ss[Oe]=L=>{S.__v_emitter=L},ss[Yl]=()=>{S.__v_emitter=void 0},ss}var Ce=typeof global<"u"?global:typeof self<"u"?self:typeof window<"u"?window:{};function $g(){return ru().__VUE_DEVTOOLS_GLOBAL_HOOK__}function ru(){return typeof navigator<"u"&&typeof window<"u"?window:typeof Ce<"u"?Ce:{}}const Bg=typeof Proxy=="function",qg="devtools-plugin:setup",zg="plugin:settings:set";let qa,Jl;function Ug(){var s;return qa!==void 0||(typeof window<"u"&&window.performance?(qa=!0,Jl=window.performance):typeof Ce<"u"&&(!((s=Ce.perf_hooks)===null||s===void 0)&&s.performance)?(qa=!0,Jl=Ce.perf_hooks.performance):qa=!1),qa}function Hg(){return Ug()?Jl.now():Date.now()}class Vg{constructor(n,a){this.target=null,this.targetQueue=[],this.onQueue=[],this.plugin=n,this.hook=a;const e={};if(n.settings)for(const o in n.settings){const p=n.settings[o];e[o]=p.defaultValue}const t=`__vue-devtools-plugin-settings__${n.id}`;let l=Object.assign({},e);try{const o=localStorage.getItem(t),p=JSON.parse(o);Object.assign(l,p)}catch{}this.fallbacks={getSettings(){return l},setSettings(o){try{localStorage.setItem(t,JSON.stringify(o))}catch{}l=o},now(){return Hg()}},a&&a.on(zg,(o,p)=>{o===this.plugin.id&&this.fallbacks.setSettings(p)}),this.proxiedOn=new Proxy({},{get:(o,p)=>this.target?this.target.on[p]:(...c)=>{this.onQueue.push({method:p,args:c})}}),this.proxiedTarget=new Proxy({},{get:(o,p)=>this.target?this.target[p]:p==="on"?this.proxiedOn:Object.keys(this.fallbacks).includes(p)?(...c)=>(this.targetQueue.push({method:p,args:c,resolve:()=>{}}),this.fallbacks[p](...c)):(...c)=>new Promise(r=>{this.targetQueue.push({method:p,args:c,resolve:r})})})}async setRealTarget(n){this.target=n;for(const a of this.onQueue)this.target.on[a.method](...a.args);for(const a of this.targetQueue)a.resolve(await this.target[a.method](...a.args))}}function Gg(s,n){const a=s,e=ru(),t=$g(),l=Bg&&a.enableEarlyProxy;if(t&&(e.__VUE_DEVTOOLS_PLUGIN_API_AVAILABLE__||!l))t.emit(qg,s,n);else{const o=l?new Vg(a,t):null;(e.__VUE_DEVTOOLS_PLUGINS__=e.__VUE_DEVTOOLS_PLUGINS__||[]).push({pluginDescriptor:a,setupFn:n,proxy:o}),o&&n(o.proxiedTarget)}}const iu="vue-i18n: composer properties",bl={"vue-devtools-plugin-vue-i18n":"Vue I18n DevTools","vue-i18n-resource-inspector":"Vue I18n DevTools","vue-i18n-timeline":"Vue I18n"},Wg={"vue-i18n-resource-inspector":"Search for scopes ..."},Kg={"vue-i18n-timeline":16764185};let Zl;async function Xg(s,n){return new Promise((a,e)=>{try{Gg({id:"vue-devtools-plugin-vue-i18n",label:bl["vue-devtools-plugin-vue-i18n"],packageName:"vue-i18n",homepage:"https://vue-i18n.intlify.dev",logo:"https://vue-i18n.intlify.dev/vue-i18n-devtools-logo.png",componentStateTypes:[iu],app:s},t=>{Zl=t,t.on.visitComponentTree(({componentInstance:o,treeNode:p})=>{Yg(o,p,n)}),t.on.inspectComponent(({componentInstance:o,instanceData:p})=>{o.vnode.el&&o.vnode.el.__VUE_I18N__&&p&&Qg(p,o.vnode.el.__VUE_I18N__)}),t.addInspector({id:"vue-i18n-resource-inspector",label:bl["vue-i18n-resource-inspector"],icon:"language",treeFilterPlaceholder:Wg["vue-i18n-resource-inspector"]}),t.on.getInspectorTree(o=>{o.app===s&&o.inspectorId==="vue-i18n-resource-inspector"&&ab(o,n)});const l=new Map;t.on.getInspectorState(async o=>{if(o.app===s&&o.inspectorId==="vue-i18n-resource-inspector")if(t.unhighlightElement(),tb(o,n),o.nodeId==="global"){if(!l.has(o.app)){const[p]=await t.getComponentInstances(o.app);l.set(o.app,p)}t.highlightElement(l.get(o.app))}else{const p=eb(o.nodeId,n);p&&t.highlightElement(p)}}),t.on.editInspectorState(o=>{o.app===s&&o.inspectorId==="vue-i18n-resource-inspector"&&ob(o,n)}),t.addTimelineLayer({id:"vue-i18n-timeline",label:bl["vue-i18n-timeline"],color:Kg["vue-i18n-timeline"]}),a(!0)})}catch(t){console.error(t),e(!1)}})}function uu(s){return s.type.name||s.type.displayName||s.type.__file||"Anonymous"}function Yg(s,n,a){const e=a.global;if(s&&s.vnode.el&&s.vnode.el.__VUE_I18N__&&s.vnode.el.__VUE_I18N__!==e){const t={label:`i18n (${uu(s)} Scope)`,textColor:0,backgroundColor:16764185};n.tags.push(t)}}function Qg(s,n){const a=iu;s.state.push({type:a,key:"locale",editable:!0,value:n.locale.value}),s.state.push({type:a,key:"availableLocales",editable:!1,value:n.availableLocales}),s.state.push({type:a,key:"fallbackLocale",editable:!0,value:n.fallbackLocale.value}),s.state.push({type:a,key:"inheritLocale",editable:!0,value:n.inheritLocale}),s.state.push({type:a,key:"messages",editable:!1,value:Fo(n.messages.value)}),s.state.push({type:a,key:"datetimeFormats",editable:!1,value:n.datetimeFormats.value}),s.state.push({type:a,key:"numberFormats",editable:!1,value:n.numberFormats.value})}function Fo(s){const n={};return Object.keys(s).forEach(a=>{const e=s[a];Os(e)&&"source"in e?n[a]=nb(e):En(e)&&e.loc&&e.loc.source?n[a]=e.loc.source:vs(e)?n[a]=Fo(e):n[a]=e}),n}const Jg={"<":"&lt;",">":"&gt;",'"':"&quot;","&":"&amp;"};function Zg(s){return s.replace(/[<>"&]/g,sb)}function sb(s){return Jg[s]||s}function nb(s){return{_custom:{type:"function",display:`<span>ƒ</span> ${s.source?`("${Zg(s.source)}")`:"(?)"}`}}}function ab(s,n){s.rootNodes.push({id:"global",label:"Global Scope"});const a=n.global;for(const[e,t]of n.__instances){const l=t;a!==l&&s.rootNodes.push({id:l.id.toString(),label:`${uu(e)} Scope`})}}function eb(s,n){let a=null;if(s!=="global"){for(const[e,t]of n.__instances.entries())if(t.id.toString()===s){a=e;break}}return a}function hu(s,n){if(s==="global")return n.global;{const a=Array.from(n.__instances.values()).find(e=>e.id.toString()===s);return a||null}}function tb(s,n){const a=hu(s.nodeId,n);return a&&(s.state=lb(a)),null}function lb(s){const n={},a="Locale related info",e=[{type:a,key:"locale",editable:!0,value:s.locale.value},{type:a,key:"fallbackLocale",editable:!0,value:s.fallbackLocale.value},{type:a,key:"availableLocales",editable:!1,value:s.availableLocales},{type:a,key:"inheritLocale",editable:!0,value:s.inheritLocale}];n[a]=e;const t="Locale messages info",l=[{type:t,key:"messages",editable:!1,value:Fo(s.messages.value)}];n[t]=l;{const o="Datetime formats info",p=[{type:o,key:"datetimeFormats",editable:!1,value:s.datetimeFormats.value}];n[o]=p;const c="Datetime formats info",r=[{type:c,key:"numberFormats",editable:!1,value:s.numberFormats.value}];n[c]=r}return n}function so(s,n){if(Zl){let a;n&&"groupId"in n&&(a=n.groupId,delete n.groupId),Zl.addTimelineEvent({layerId:"vue-i18n-timeline",event:{title:s,groupId:a,time:Date.now(),meta:{},data:n||{},logType:s==="compile-error"?"error":s==="fallback"||s==="missing"?"warning":"default"}})}}function ob(s,n){const a=hu(s.nodeId,n);if(a){const[e]=s.path;e==="locale"&&ns(s.state.value)?a.locale.value=s.state.value:e==="fallbackLocale"&&(ns(s.state.value)||Gs(s.state.value)||vs(s.state.value))?a.fallbackLocale.value=s.state.value:e==="inheritLocale"&&Ns(s.state.value)&&(a.inheritLocale=s.state.value)}}const $o={tag:{type:[String,Object]},locale:{type:String},scope:{type:String,validator:s=>s==="parent"||s==="global",default:"parent"},i18n:{type:Object}};function pb({slots:s},n){return n.length===1&&n[0]==="default"?(s.default?s.default():[]).reduce((e,t)=>[...e,...t.type===js?t.children:[t]],[]):n.reduce((a,e)=>{const t=s[e];return t&&(a[e]=t()),a},Rs())}function du(){return js}function cb(s){return Gs(s)&&!ns(s[0])}function mu(s,n,a,e){const{slots:t,attrs:l}=n;return()=>{const o={part:!0};let p=Rs();s.locale&&(o.locale=s.locale),ns(s.format)?o.key=s.format:vs(s.format)&&(ns(s.format.key)&&(o.key=s.format.key),p=Object.keys(s.format).reduce((h,d)=>a.includes(d)?en(Rs(),h,{[d]:s.format[d]}):h,Rs()));const c=e(s.value,o,p);let r=[o.key];Gs(c)?r=c.map((h,d)=>{const f=t[h.type],m=f?f({[h.type]:h.value,index:d,parts:c}):[h.value];return cb(m)&&(m[0].key=`${h.type}-${d}`),m}):ns(c)&&(r=[c]);const i=en(Rs(),l),u=ns(s.tag)||vs(s.tag)?s.tag:du();return Ve(u,i,r)}}const rb=Es({name:"i18n-d",props:en({value:{type:[Number,Date],required:!0},format:{type:[String,Object]}},$o),setup(s,n){const a=s.i18n||Hs({useScope:s.scope,__useComponent:!0});return mu(s,n,au,(...e)=>a[Kl](...e))}}),Lc=rb,ib=Es({name:"i18n-n",props:en({value:{type:Number,required:!0},format:{type:[String,Object]}},$o),setup(s,n){const a=s.i18n||Hs({useScope:s.scope,__useComponent:!0});return mu(s,n,eu,(...e)=>a[Xl](...e))}}),Mc=ib,ub=Es({name:"i18n-t",props:en({},{keypath:{type:String,required:!0},plural:{type:[Number,String],validator:s=>Xs(s)||!isNaN(s)}},$o),setup(s,n){const{slots:a,attrs:e}=n,t=s.i18n||Hs({useScope:s.scope,__useComponent:!0});return()=>{const l=Object.keys(a).filter(u=>u[0]!=="_"),o=Rs();s.locale&&(o.locale=s.locale),s.plural!==void 0&&(o.plural=ns(s.plural)?+s.plural:s.plural);const p=pb(n,l),c=t[Wl](s.keypath,p,o),r=en(Rs(),e),i=ns(s.tag)||vs(s.tag)?s.tag:du();return Ve(i,r,c)}}}),Dc=ub;function hb(s,...n){const a=ys(n[0])?n[0]:{};(Ns(a.globalInstall)?a.globalInstall:!0)&&([Dc.name,"I18nT"].forEach(t=>s.component(t,Dc)),[Mc.name,"I18nN"].forEach(t=>s.component(t,Mc)),[Lc.name,"I18nD"].forEach(t=>s.component(t,Lc)))}const db=Un("global-vue-i18n");function mb(s={}){const n=Ns(s.globalInjection)?s.globalInjection:!0,a=new Map,[e,t]=jb(s),l=Un("vue-i18n");function o(i){return a.get(i)||null}function p(i,u){a.set(i,u)}function c(i){a.delete(i)}const r={async install(i,...u){if(i.__VUE_I18N__=r,i.__VUE_I18N_SYMBOL__=l,i.provide(i.__VUE_I18N_SYMBOL__,r),ys(u[0])){const f=u[0];r.__composerExtend=f.__composerExtend}let h=null;n&&(h=Cb(i,r.global)),hb(i,...u);const d=i.unmount;i.unmount=()=>{h&&h(),r.dispose(),d()};{if(!await Xg(i,r))throw zn(Ks.CANNOT_SETUP_VUE_DEVTOOLS_PLUGIN);const m=Ni(),k=t;k[Oe]&&k[Oe](m),m.on("*",so)}},get global(){return t},dispose(){e.stop()},__instances:a,__getInstance:o,__setInstance:p,__deleteInstance:c};return r}function Hs(s={}){const n=ca();if(n==null)throw zn(Ks.MUST_BE_CALL_SETUP_TOP);if(!n.isCE&&n.appContext.app!=null&&!n.appContext.app.__VUE_I18N_SYMBOL__)throw zn(Ks.NOT_INSTALLED);const a=fb(n),e=bb(a),t=pu(n),l=gb(s,t);if(l==="global")return Ig(e,s,t),e;if(l==="parent"){let c=_b(a,n,s.__useComponent);return c==null&&(ka(No(Qa.NOT_FOUND_PARENT_SCOPE)),c=e),c}const o=a;let p=o.__getInstance(n);if(p==null){const c=en({},s);"__i18n"in t&&(c.__i18n=t.__i18n),e&&(c.__root=e),p=cu(c),o.__composerExtend&&(p[Ql]=o.__composerExtend(p)),vb(o,n,p),o.__setInstance(n,p)}else if(l==="local")throw zn(Ks.DUPLICATE_USE_I18N_CALLING);return p}function jb(s){const n=ro(),a=n.run(()=>cu(s));if(a==null)throw zn(Ks.UNEXPECTED_ERROR);return[n,a]}function fb(s){const n=_n(s.isCE?db:s.appContext.app.__VUE_I18N_SYMBOL__);if(!n)throw zn(s.isCE?Ks.NOT_INSTALLED_WITH_PROVIDE:Ks.UNEXPECTED_ERROR);return n}function gb(s,n){return Vt(s)?"__i18n"in n?"local":"global":s.useScope?s.useScope:"local"}function bb(s){return s.global}function _b(s,n,a=!1){let e=null;const t=n.root;let l=yb(n,a);for(;l!=null&&(e=s.__getInstance(l),!(e!=null||t===l));)l=l.parent;return e}function yb(s,n=!1){return s==null?null:n&&s.vnode.ctx||s.parent}function vb(s,n,a){let e=null;kn(()=>{if(n.vnode.el){n.vnode.el.__VUE_I18N__=a,e=Ni();const t=a;t[Oe]&&t[Oe](e),e.on("*",so)}},n),Co(()=>{const t=a;n.vnode.el&&n.vnode.el.__VUE_I18N__&&(e&&e.off("*",so),t[Yl]&&t[Yl](),delete n.vnode.el.__VUE_I18N__),s.__deleteInstance(n);const l=t[Ql];l&&(l(),delete t[Ql])},n)}const wb=["locale","fallbackLocale","availableLocales"],Oc=["t","rt","d","n","tm","te"];function Cb(s,n){const a=Object.create(null);return wb.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l)throw zn(Ks.UNEXPECTED_ERROR);const o=zs(l.value)?{get(){return l.value.value},set(p){l.value.value=p}}:{get(){return l.get&&l.get()}};Object.defineProperty(a,t,o)}),s.config.globalProperties.$i18n=a,Oc.forEach(t=>{const l=Object.getOwnPropertyDescriptor(n,t);if(!l||!l.value)throw zn(Ks.UNEXPECTED_ERROR);Object.defineProperty(s.config.globalProperties,`$${t}`,l)}),()=>{delete s.config.globalProperties.$i18n,Oc.forEach(t=>{delete s.config.globalProperties[`$${t}`]})}}rg(qf);ig(tg);ug(Wi);{const s=tf();s.__INTLIFY__=!0,zf(s.__INTLIFY_DEVTOOLS_GLOBAL_HOOK__)}Rg();const kb={tags:"标签",articles:"文章",words:"字数",prevPage:"上一页",nextPage:"下一页",pagination:"分页导航",search:"搜索",theme:"主题",language:"语言",menu:"菜单",close:"关闭",searchPlaceholder:"搜索文章标题或正文...",searchNoResults:"未找到匹配的文章",searchUnavailable:"搜索索引加载失败",searchEscHint:"Esc 关闭",searchResultsLabel:"条结果",categories:"分类",resources:"资源",about:"关于",greeting:"你好，我是",greetingPrefix:"//",developer:"开发者",wordUnit:"字",recentPosts:"最近文章",noPosts:"暂无文章",notes:"笔记",projects:"项目",topics:"课题",seeMore:"查看更多",countPosts:"{count} 篇",countProjects:"{count} 个项目",countTopics:"{count} 个课题",countWords:"{count} 字",copyCode:"复制代码",copyFailed:"复制失败",copyArticle:"复制文章",copyTable:"复制表格",copied:"已复制",anchorHeading:"置顶当前标题",uncategorized:"未分类",updatedAt:"更新于",postReadingTime:"{minutes} 分钟",articleReadingTime:"阅读约 {minutes} 分钟",tableOfContents:"本页目录",openToc:"打开目录",backToTop:"回到顶部",resourceSubtitle:"生物信息学与结构生物学领域常用工具",pageNotFound:"页面不存在或已被移动",backHome:"返回首页",metaHome:"首页",metaCategories:"分类",metaResources:"资源",metaAbout:"关于",metaNotFound:"404",experience:"经历",thanks:"感谢您的关注！",designedByPrefix:"由",designedBySuffix:"设计",clearFilter:"清除筛选",backToArticle:"返回文章"},xb={tags:"Tags",articles:"Articles",words:"Words",prevPage:"Previous",nextPage:"Next",pagination:"Pagination",search:"Search",theme:"Theme",language:"Language",menu:"Menu",close:"Close",searchPlaceholder:"Search articles by title or content...",searchNoResults:"No matching articles found",searchUnavailable:"Search index failed to load",searchEscHint:"Esc to close",searchResultsLabel:"results",categories:"Categories",resources:"Resources",about:"About",greeting:"Hello, I'm",greetingPrefix:"//",developer:"Developer",wordUnit:"words",recentPosts:"Recent Posts",noPosts:"No posts yet",notes:"Notes",projects:"Projects",topics:"Topics",seeMore:"See More",countPosts:"{count} posts",countProjects:"{count} projects",countTopics:"{count} topics",countWords:"{count} words",copyCode:"Copy code",copyFailed:"Copy failed",copyArticle:"Copy article",copyTable:"Copy table",copied:"Copied",anchorHeading:"Anchor to heading",uncategorized:"Uncategorized",updatedAt:"Updated at",postReadingTime:"{minutes} min",articleReadingTime:"Reading about {minutes} minutes",tableOfContents:"On this page",openToc:"Open table of contents",backToTop:"Back to top",resourceSubtitle:"Commonly used tools in bioinformatics and structural biology",pageNotFound:"This page could not be found",backHome:"Back to home",metaHome:"Home",metaCategories:"Categories",metaResources:"Resources",metaAbout:"About",metaNotFound:"404",experience:"Experience",thanks:"Thank you for your attention!",designedByPrefix:"Designed by",designedBySuffix:"",clearFilter:"Clear filter",backToArticle:"Back to article"},jn={url:"https://zorrooz.github.io",author:"zorrooz",email:"zorrooz@163.com",github:"https://github.com/zorrooz",startYear:2025},Ic="zorrooz’s blog",Zn=["zh-CN","en-US"],Pb=["zh","en"],_l=["auto","light","dark"],ju="zh-CN",Bo="locale",no="theme",Nc="-en",ta="/article",qo=88,Sb={zh:"zh-CN",en:"en-US"},ao={"zh-CN":"zh","en-US":"en"},Et={"zh-CN":"zh-CN","en-US":"en"},Kt=s=>ao[s],zo=s=>{const n=s.match(/^\/(zh|en)(?=\/|$)/);return n?Sb[n[1]]:null},fu=s=>s.replace(/^\/(zh|en)(?=\/|$)/,""),Na=s=>s===Zn[1]?Zn[1]:Zn[0],za=()=>{if(typeof window<"u"){const s=localStorage.getItem(Bo);if(s===Zn[1])return ao[Zn[1]];if(s===Zn[0])return ao[Zn[0]];if(navigator.language&&navigator.language.toLowerCase().startsWith("en"))return"en"}return"zh"},Eb=()=>{const s=globalThis.__GBLOG_LOCALE__;return s===Zn[0]||s===Zn[1]?s:null},ee=()=>typeof window<"u"?Eb()??Na(localStorage.getItem(Bo)):ju,Tb=s=>{typeof window<"u"&&localStorage.setItem(Bo,s)},gu=ee();typeof document<"u"&&(document.documentElement.lang=Et[gu]);const Ne=mb({locale:gu,fallbackLocale:ju,messages:{"zh-CN":kb,"en-US":xb},globalInjection:!0}),Uo=ji("locale",()=>{const s=ps(ee());return{locale:s,setLocale:e=>{s.value=e,Ne.global.locale.value=e,Tb(e),typeof document<"u"&&(document.documentElement.lang=Et[e])},initLocale:()=>{if(typeof window<"u"){const e=zo(window.location.pathname);e&&(s.value=e)}Ne.global.locale.value=s.value,typeof document<"u"&&(document.documentElement.lang=Et[s.value])}}}),Fc=s=>{if(typeof window>"u"||typeof document>"u")return;const n=document.documentElement;if(s==="auto"){const a=window.matchMedia("(prefers-color-scheme: dark)").matches;n.setAttribute("data-bs-theme",a?"dark":"light"),localStorage.removeItem(no)}else n.setAttribute("data-bs-theme",s),localStorage.setItem(no,s)},Ab=()=>{if(typeof window>"u")return"auto";const s=localStorage.getItem(no);return s==="light"||s==="dark"?s:"auto"},Ho=ji("theme",()=>{const s=ps(Ab());return{theme:s,toggleTheme:()=>{const t=(_l.indexOf(s.value)+1)%_l.length;s.value=_l[t],Fc(s.value)},initTheme:()=>{Fc(s.value)}}}),bu=()=>typeof window<"u"&&typeof window.matchMedia=="function"&&window.matchMedia("(prefers-reduced-motion: reduce)").matches,Fa=(s="smooth")=>{typeof window>"u"||window.scrollTo({top:0,behavior:bu()?"auto":s})},Rb=()=>typeof document<"u"&&!!document.body?.style&&document.body.style.overflow==="hidden",_u=(s,n=0,a={})=>{const e=Math.max(0,s.getBoundingClientRect().top+window.scrollY-n),t=()=>{window.scrollTo({top:e,behavior:bu()?"auto":"smooth"})};try{Rb()?setTimeout(t,a.lockedDelay??80):t()}catch{t()}};function ne(s,n){const a=document.documentElement;a&&(a.style.overflow=s,a.style.overscrollBehavior=n)}function Xt(s,n){const a=document.body;a&&(a.style.overflow=s,a.style.overscrollBehavior=n)}function Lb(){ne("hidden","contain"),Xt("hidden","contain")}function $c(){ne("",""),Xt("","")}function Mb(){const s=window.scrollY||window.pageYOffset||document.documentElement.scrollTop||0;try{const n=document.body;n&&(n.style.position="fixed",n.style.top=`-${s}px`,n.style.left="0",n.style.right="0",n.style.overflow="hidden"),ne("","contain")}catch{ne("hidden","contain"),Xt("hidden","contain")}return s}function Bc(s){try{const n=document.body;n&&(n.style.position="",n.style.top="",n.style.left="",n.style.right="",n.style.overflow=""),ne("",""),typeof s=="number"&&window.scrollTo(0,s)}catch{ne("",""),Xt("","")}}const eo=s=>{const n=s.replace(/^\/+/,"").split("/"),a=n[0]===ta.slice(1)?1:0;return`${n[a]}/${n[a+1]}/${n.slice(a+2).join("/")}.md`},Db=s=>`${ta}/${s.replace(/\.md$/,"")}`,Fe=s=>ga(Db(s)),Tt=s=>s.replace(/\.md$/,"").replace(/-en$/,""),dt=s=>Array.isArray(s)?s.join("/"):typeof s=="string"&&s.length?s:"",ga=s=>`/${Kt(Na(Ne.global.locale.value))}${s.startsWith("/")?s:`/${s}`}`,Ob=(s,n)=>{const e=(zo(n.path)??Na(Ne.global.locale.value))==="zh-CN"?"en-US":"zh-CN";s.push({path:`/${Kt(e)}${fu(n.path)}`,query:n.query})};function qc(s,n,a){n&&(s.push({path:`/${Kt(a.locale)}/`,query:{...a.query??{},tag:n}}).catch(()=>{}),a.scroll!==!1&&vn(()=>Fa()))}function wa(s,n){const a=ps(n),{locale:e}=Hs(),t=()=>{try{a.value=s()}catch(l){console.error("Failed to load localized content:",l),a.value=n}};return t(),Ws(e,t),{data:a,reload:t}}const yu=`Hello, I am zorrooz, currently focused on structural analysis and functional studies of biological macromolecules, while also exploring computational biology.
`,vu=[{year:"2021 – 2025",title:"Lanzhou University · B.S. in Bioinformatics",desc:"Basic programming training and foundational knowledge in biology"},{year:"2025 – present",title:"Lanzhou University · M.S. in Biology",desc:"Structural and functional studies of biological macromolecules"}],wu=[{title:"Common Tools",items:[{name:"Python · R · Bash · Git",desc:"Primary languages and version control tools for daily research and development"},{name:"Nextflow · Snakemake",desc:"Bioinformatics pipeline orchestration and workflow management"},{name:"RELION · cryoSPARC",desc:"Cryo-EM single-particle reconstruction and 3D structure validation"},{name:"VS Code · Ubuntu · WSL",desc:"Development environment and terminal workflow"}]},{title:"Fields",items:[{name:"Structural Biology",desc:"Cryo-EM"},{name:"Computational Biology",desc:"Binder design and virtual screening"},{name:"Programming",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],Cu=[{label:"Email",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],Ib={introduction:yu,experience:vu,section:wu,contacts:Cu},Nb=Object.freeze(Object.defineProperty({__proto__:null,contacts:Cu,default:Ib,experience:vu,introduction:yu,section:wu},Symbol.toStringTag,{value:"Module"})),ku=`你好，我是 zorrooz，当前专注于生物大分子的结构解析与功能研究，同时涉猎计算生物学。
`,xu=[{year:"2021 – 2025",title:"兰州大学 本科 · 生物信息学",desc:"基本编程训练与生物学基础知识"},{year:"2025 – 至今",title:"兰州大学 硕士 · 生物学",desc:"生物大分子的结构与功能研究"}],Pu=[{title:"常用工具",items:[{name:"Python · R · Bash · Git",desc:"日常研究与开发的主力语言与版本控制工具"},{name:"Nextflow · Snakemake",desc:"生信流水线编排与工作流管理"},{name:"RELION · cryoSPARC",desc:"冷冻电镜单颗粒重构与三维结构验证"},{name:"VS Code · Ubuntu · WSL",desc:"开发环境与终端工作流"}]},{title:"领域",items:[{name:"结构生物学",desc:"冷冻电镜"},{name:"计算生物学",desc:"binder 设计与虚拟筛选"},{name:"编程",desc:"R · Python · C · Perl · Bash · JavaScript"}]}],Su=[{label:"邮箱",value:"zorrooz@163.com",link:"mailto:zorrooz@163.com",icon:"fas fa-envelope"},{label:"GitHub",value:"zorrooz",link:"https://github.com/zorrooz",icon:"fab fa-github"}],Fb={introduction:ku,experience:xu,section:Pu,contacts:Su},$b=Object.freeze(Object.defineProperty({__proto__:null,contacts:Su,default:Fb,experience:xu,introduction:ku,section:Pu},Symbol.toStringTag,{value:"Module"})),Bb=JSON.parse('[{"title":"Notes","items":[{"name":"ComputationalBiology","title":"Computational Biology","desc":"Mainstream tools for protein design and virtual screening","tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology","Virtual Screening","Molecular Docking","AutoDock Vina"],"stats":{"postsCount":2,"totalWords":1211,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"Protein Design","articles":[{"title":"Mainstream Tools for Protein Design","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools-en","wordCount":618,"tags":["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":618,"latestDate":"2026-08-04"}},{"key":"virtual-screening","title":"Virtual Screening","articles":[{"title":"Mainstream Tools for Virtual Screening","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools-en","wordCount":593,"tags":["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"]}],"stats":{"postsCount":1,"totalWords":593,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools-en"},{"name":"Programming","title":"Programming","desc":"Detailed tutorials on R, Python, Linux, and Bash programming","tags":["Bash","Shell","Scripting","Tutorial","Linux","Command Line","Python","Advanced","Getting Started","NumPy","Pandas","Data Processing","R Language","ggplot2","Data Visualization","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":1972,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python Advanced: Functions, Classes, and Modules","articleUrl":"/article/notes/Programming/python/python-advanced-en","wordCount":244,"tags":["Python","Advanced","Tutorial"]},{"title":"Introduction to Python Programming: Environment, Syntax, and Data Types","articleUrl":"/article/notes/Programming/python/python-basics-en","wordCount":333,"tags":["Python","Getting Started","Tutorial"]},{"title":"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data-en","wordCount":314,"tags":["Python","NumPy","Pandas","Data Processing","Tutorial"]}],"stats":{"postsCount":3,"totalWords":891,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R Language Primer: Data Structures, Vectorization, and Functions","articleUrl":"/article/notes/Programming/r/r-basics-en","wordCount":256,"tags":["R Language","Getting Started","Tutorial"]},{"title":"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts","articleUrl":"/article/notes/Programming/r/r-ggplot2-en","wordCount":157,"tags":["R Language","ggplot2","Data Visualization","Tutorial"]},{"title":"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning","articleUrl":"/article/notes/Programming/r/r-tidyverse-en","wordCount":245,"tags":["R Language","tidyverse","dplyr","Tutorial"]}],"stats":{"postsCount":3,"totalWords":658,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux Command Line Basics: File System, Permissions, and Text Processing","articleUrl":"/article/notes/Programming/linux/linux-basics-en","wordCount":238,"tags":["Linux","Command Line","Tutorial"]}],"stats":{"postsCount":1,"totalWords":238,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts","articleUrl":"/article/notes/Programming/bash/bash-scripting-en","wordCount":185,"tags":["Bash","Shell","Scripting","Tutorial"]}],"stats":{"postsCount":1,"totalWords":185,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced-en"},{"name":"StructuralBiology","title":"Structural Biology","desc":"Cryo-EM data processing, visualization of biomacromolecular structures, and atomic modeling","tags":["Cryo-Electron Microscopy","cryo-EM","Review","Data Processing","RELION","Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial","Structure Visualization","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":2670,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"Cryo-EM","articles":[{"title":"Review of Cryo-EM Technology","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview-en","wordCount":977,"tags":["Cryo-Electron Microscopy","cryo-EM","Review"]},{"title":"Cryo-EM Single Particle Analysis: Full Data Processing Workflow","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow-en","wordCount":765,"tags":["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"]}],"stats":{"postsCount":2,"totalWords":1742,"latestDate":"2026-08-04"}},{"key":"visualization","title":"Structural Visualization","articles":[{"title":"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization-en","wordCount":362,"tags":["Structure Visualization","PyMOL","ChimeraX","Tutorial"]}],"stats":{"postsCount":1,"totalWords":362,"latestDate":"2026-08-04"}},{"key":"modeling","title":"Atomic Modeling","articles":[{"title":"Atomic Modeling and Refinement","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling-en","wordCount":566,"tags":["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"]}],"stats":{"postsCount":1,"totalWords":566,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview-en"}]},{"title":"Projects","items":[{"name":"Plotit","title":"plotit","desc":"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage","tags":["R","ggplot2","Data Visualization","Declarative API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management","tags":["Vue","Vite","Blog","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"Topics","items":[{"name":"StructureOfProteinDemo","title":"Protein Structure Determination Demo Project","desc":"Protein structure placeholder demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),qb=Object.freeze(Object.defineProperty({__proto__:null,default:Bb},Symbol.toStringTag,{value:"Module"})),zb=JSON.parse('[{"title":"笔记","items":[{"name":"ComputationalBiology","title":"计算生物学","desc":"蛋白质设计与虚拟筛选的主流工具","tags":["虚拟筛选","分子对接","AutoDock Vina","计算生物学","蛋白质设计","Rosetta","AlphaFold","RFdiffusion"],"stats":{"postsCount":2,"totalWords":1899,"latestDate":"2026-08-04"},"categories":[{"key":"protein-design","title":"蛋白质设计","articles":[{"title":"蛋白质设计主流工具","articleUrl":"/article/notes/ComputationalBiology/protein-design/protein-design-tools","wordCount":1002,"tags":["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"]}],"stats":{"postsCount":1,"totalWords":1002,"latestDate":""}},{"key":"virtual-screening","title":"虚拟筛选","articles":[{"title":"虚拟筛选主流工具","articleUrl":"/article/notes/ComputationalBiology/virtual-screening/virtual-screening-tools","wordCount":897,"tags":["虚拟筛选","分子对接","AutoDock Vina","计算生物学"]}],"stats":{"postsCount":1,"totalWords":897,"latestDate":"2026-08-04"}}],"root":"/article/notes/ComputationalBiology/protein-design/protein-design-tools"},{"name":"Programming","title":"编程","desc":"R、Python、Linux 与 Bash 编程的详细教程","tags":["Bash","Shell","脚本编程","教程","Linux","命令行","Python","进阶","入门","NumPy","Pandas","数据处理","R语言","ggplot2","数据可视化","tidyverse","dplyr"],"stats":{"postsCount":8,"totalWords":2903,"latestDate":"2026-08-04"},"categories":[{"key":"python","title":"Python","articles":[{"title":"Python 进阶：函数、类与模块","articleUrl":"/article/notes/Programming/python/python-advanced","wordCount":376,"tags":["Python","进阶","教程"]},{"title":"Python 编程入门：环境、语法与数据类型","articleUrl":"/article/notes/Programming/python/python-basics","wordCount":495,"tags":["Python","入门","教程"]},{"title":"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas","articleUrl":"/article/notes/Programming/python/python-data","wordCount":473,"tags":["Python","NumPy","Pandas","数据处理","教程"]}],"stats":{"postsCount":3,"totalWords":1344,"latestDate":"2026-08-04"}},{"key":"r","title":"R","articles":[{"title":"R 语言入门：数据结构、向量化与函数","articleUrl":"/article/notes/Programming/r/r-basics","wordCount":403,"tags":["R语言","入门","教程"]},{"title":"R ggplot2 数据可视化：图层语法与常用图表","articleUrl":"/article/notes/Programming/r/r-ggplot2","wordCount":245,"tags":["R语言","ggplot2","数据可视化","教程"]},{"title":"R tidyverse 数据操作：dplyr 管道与数据清洗","articleUrl":"/article/notes/Programming/r/r-tidyverse","wordCount":324,"tags":["R语言","tidyverse","dplyr","教程"]}],"stats":{"postsCount":3,"totalWords":972,"latestDate":"2026-08-04"}},{"key":"linux","title":"Linux","articles":[{"title":"Linux 命令行基础：文件系统、权限与文本处理","articleUrl":"/article/notes/Programming/linux/linux-basics","wordCount":328,"tags":["Linux","命令行","教程"]}],"stats":{"postsCount":1,"totalWords":328,"latestDate":"2026-08-04"}},{"key":"bash","title":"Bash","articles":[{"title":"Bash 编程：变量、条件、循环与实用脚本","articleUrl":"/article/notes/Programming/bash/bash-scripting","wordCount":259,"tags":["Bash","Shell","脚本编程","教程"]}],"stats":{"postsCount":1,"totalWords":259,"latestDate":"2026-08-04"}}],"root":"/article/notes/Programming/python/python-advanced"},{"name":"StructuralBiology","title":"结构生物学","desc":"冷冻电镜数据处理、生物大分子结构可视化与原子建模","tags":["冷冻电镜","cryo-EM","综述","数据处理","RELION","原子建模","Coot","phenix","结构精修","教程","结构可视化","PyMOL","ChimeraX"],"stats":{"postsCount":4,"totalWords":3960,"latestDate":"2026-08-04"},"categories":[{"key":"cryoem","title":"冷冻电镜","articles":[{"title":"冷冻电镜技术综述","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-overview","wordCount":1485,"tags":["冷冻电镜","cryo-EM","综述"]},{"title":"冷冻电镜单颗粒分析：数据处理全流程","articleUrl":"/article/notes/StructuralBiology/cryoem/cryoem-workflow","wordCount":1111,"tags":["冷冻电镜","cryo-EM","数据处理","RELION"]}],"stats":{"postsCount":2,"totalWords":2596,"latestDate":"2026-08-04"}},{"key":"visualization","title":"结构可视化","articles":[{"title":"生物大分子结构可视化：PyMOL 与 ChimeraX 实战","articleUrl":"/article/notes/StructuralBiology/visualization/structure-visualization","wordCount":529,"tags":["结构可视化","PyMOL","ChimeraX","教程"]}],"stats":{"postsCount":1,"totalWords":529,"latestDate":"2026-08-04"}},{"key":"modeling","title":"原子建模","articles":[{"title":"原子建模与精修","articleUrl":"/article/notes/StructuralBiology/modeling/atomic-modeling","wordCount":835,"tags":["原子建模","Coot","phenix","结构精修","教程"]}],"stats":{"postsCount":1,"totalWords":835,"latestDate":"2026-08-04"}}],"root":"/article/notes/StructuralBiology/cryoem/cryoem-overview"}]},{"title":"项目","items":[{"name":"Plotit","title":"plotit","desc":"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段","tags":["R","ggplot2","数据可视化","声明式API"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-12-30"},"categories":[],"root":"","github":"https://github.com/zorrooz/plotit","status":"active","language":"R","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/plotit"},{"name":"ZorroozBlog","title":"zorrooz.github.io","desc":"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理","tags":["Vue","Vite","博客","Markdown"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2025-08-29"},"categories":[],"root":"","github":"https://github.com/zorrooz/zorrooz.github.io","status":"active","language":"Vue","stars":0,"license":"MIT","version":"","url":"https://zorrooz.github.io/"}]},{"title":"课题","items":[{"name":"StructureOfProteinDemo","title":"蛋白质结构解析演示课题","desc":"蛋白质结构占位demo","tags":["Protein","Structure","Cryo-EM"],"stats":{"postsCount":0,"totalWords":0,"latestDate":"2026-04-12"},"categories":[],"root":"","doi":"10.1234/demo-structure.2026","status":"completed","journal":"Journal of Demo Science","year":2026,"authors":["Demo Author A","Demo Author B"],"url":"https://www.demo-structure.org"}]}]'),Ub=Object.freeze(Object.defineProperty({__proto__:null,default:zb},Symbol.toStringTag,{value:"Module"})),Hb=[{title:"Mainstream Tools for Protein Design",date:"2026-08-04",tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],draft:!1,description:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",relativePath:"ComputationalBiology/protein-design/protein-design-tools-en",wordCount:618,tagCount:5},{title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],draft:!1,description:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools-en",wordCount:593,tagCount:4},{title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",tags:["Bash","Shell","Scripting","Tutorial"],draft:!1,description:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",relativePath:"Programming/bash/bash-scripting-en",wordCount:185,tagCount:4},{title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",tags:["Linux","Command Line","Tutorial"],draft:!1,description:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",relativePath:"Programming/linux/linux-basics-en",wordCount:238,tagCount:3},{title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",tags:["Python","Advanced","Tutorial"],draft:!1,description:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",relativePath:"Programming/python/python-advanced-en",wordCount:244,tagCount:3},{title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",tags:["Python","Getting Started","Tutorial"],draft:!1,description:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",relativePath:"Programming/python/python-basics-en",wordCount:333,tagCount:3},{title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],draft:!1,description:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",relativePath:"Programming/python/python-data-en",wordCount:314,tagCount:5},{title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",tags:["R Language","Getting Started","Tutorial"],draft:!1,description:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",relativePath:"Programming/r/r-basics-en",wordCount:256,tagCount:3},{title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",tags:["R Language","ggplot2","Data Visualization","Tutorial"],draft:!1,description:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",relativePath:"Programming/r/r-ggplot2-en",wordCount:157,tagCount:4},{title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",tags:["R Language","tidyverse","dplyr","Tutorial"],draft:!1,description:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",relativePath:"Programming/r/r-tidyverse-en",wordCount:245,tagCount:4},{title:"Review of Cryo-EM Technology",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Review"],draft:!1,description:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",relativePath:"StructuralBiology/cryoem/cryoem-overview-en",wordCount:977,tagCount:3},{title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],draft:!1,description:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",relativePath:"StructuralBiology/cryoem/cryoem-workflow-en",wordCount:765,tagCount:4},{title:"Atomic Modeling and Refinement",date:"2026-08-04",tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],draft:!1,description:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",relativePath:"StructuralBiology/modeling/atomic-modeling-en",wordCount:566,tagCount:5},{title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],draft:!1,description:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",relativePath:"StructuralBiology/visualization/structure-visualization-en",wordCount:362,tagCount:4}],Vb=Object.freeze(Object.defineProperty({__proto__:null,default:Hb},Symbol.toStringTag,{value:"Module"})),Gb=[{title:"虚拟筛选主流工具",date:"2026-08-04",tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],draft:!1,description:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",relativePath:"ComputationalBiology/virtual-screening/virtual-screening-tools",wordCount:897,tagCount:4},{title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",tags:["Bash","Shell","脚本编程","教程"],draft:!1,description:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",relativePath:"Programming/bash/bash-scripting",wordCount:259,tagCount:4},{title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",tags:["Linux","命令行","教程"],draft:!1,description:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",relativePath:"Programming/linux/linux-basics",wordCount:328,tagCount:3},{title:"Python 进阶：函数、类与模块",date:"2026-08-04",tags:["Python","进阶","教程"],draft:!1,description:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",relativePath:"Programming/python/python-advanced",wordCount:376,tagCount:3},{title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",tags:["Python","入门","教程"],draft:!1,description:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",relativePath:"Programming/python/python-basics",wordCount:495,tagCount:3},{title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",tags:["Python","NumPy","Pandas","数据处理","教程"],draft:!1,description:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",relativePath:"Programming/python/python-data",wordCount:473,tagCount:5},{title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",tags:["R语言","入门","教程"],draft:!1,description:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",relativePath:"Programming/r/r-basics",wordCount:403,tagCount:3},{title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",tags:["R语言","ggplot2","数据可视化","教程"],draft:!1,description:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",relativePath:"Programming/r/r-ggplot2",wordCount:245,tagCount:4},{title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",tags:["R语言","tidyverse","dplyr","教程"],draft:!1,description:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",relativePath:"Programming/r/r-tidyverse",wordCount:324,tagCount:4},{title:"冷冻电镜技术综述",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","综述"],draft:!1,description:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",relativePath:"StructuralBiology/cryoem/cryoem-overview",wordCount:1485,tagCount:3},{title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",tags:["冷冻电镜","cryo-EM","数据处理","RELION"],draft:!1,description:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",relativePath:"StructuralBiology/cryoem/cryoem-workflow",wordCount:1111,tagCount:4},{title:"原子建模与精修",date:"2026-08-04",tags:["原子建模","Coot","phenix","结构精修","教程"],draft:!1,description:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",relativePath:"StructuralBiology/modeling/atomic-modeling",wordCount:835,tagCount:5},{title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",tags:["结构可视化","PyMOL","ChimeraX","教程"],draft:!1,description:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",relativePath:"StructuralBiology/visualization/structure-visualization",wordCount:529,tagCount:4},{title:"蛋白质设计主流工具",date:"",tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],draft:!1,description:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",relativePath:"ComputationalBiology/protein-design/protein-design-tools",wordCount:1002,tagCount:5}],Wb=Object.freeze(Object.defineProperty({__proto__:null,default:Gb},Symbol.toStringTag,{value:"Module"})),Kb=[{id:1,no:1,title:"Mainstream Tools for Protein Design",date:"2026-08-04",category:["Computational Biology","protein-design"],tags:["Protein Design","Rosetta","AlphaFold","RFdiffusion","Computational Biology"],preview:"A panorama of protein design tools: Rosetta energy optimization, AlphaFold structure prediction, RFdiffusion generative design, ESMFold fast prediction, and binder design workflow",wordCount:618},{id:2,no:2,title:"Mainstream Tools for Virtual Screening",date:"2026-08-04",category:["Computational Biology","virtual-screening"],tags:["Virtual Screening","Molecular Docking","AutoDock Vina","Computational Biology"],preview:"A comprehensive overview of virtual screening tools: docking software such as AutoDock Vina, Glide, and DOCK, receptor/ligand preparation, scoring functions, high-throughput screening workflows, and validation strategies",wordCount:593},{id:3,no:3,title:"Bash Programming: Variables, Conditionals, Loops, and Practical Scripts",date:"2026-08-04",category:["Programming","bash"],tags:["Bash","Shell","Scripting","Tutorial"],preview:"A complete Bash scripting tutorial: variables and parameters, conditional statements, loops, functions, arrays, and string processing, with practical scripts for batch processing and bioinformatics",wordCount:185},{id:4,no:4,title:"Linux Command Line Basics: File System, Permissions, and Text Processing",date:"2026-08-04",category:["Programming","linux"],tags:["Linux","Command Line","Tutorial"],preview:"Core skills for the Linux command line: file system navigation, file operations, permission management, the three musketeers of text processing (grep/sed/awk), and process management",wordCount:238},{id:5,no:5,title:"Python Advanced: Functions, Classes, and Modules",date:"2026-08-04",category:["Programming","python"],tags:["Python","Advanced","Tutorial"],preview:"Deep dive into Python functional and object-oriented programming: parameter passing, lambda, decorators, classes and inheritance, exception handling, modules and packages",wordCount:244},{id:6,no:6,title:"Introduction to Python Programming: Environment, Syntax, and Data Types",date:"2026-08-04",category:["Programming","python"],tags:["Python","Getting Started","Tutorial"],preview:"Learn Python from scratch: environment setup, basic syntax, built-in data types, and control flow, with runnable example code.",wordCount:333},{id:7,no:7,title:"Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas",date:"2026-08-04",category:["Programming","python"],tags:["Python","NumPy","Pandas","Data Processing","Tutorial"],preview:"A complete practical tutorial on file I/O, regular expressions, and the scientific computing stack (NumPy/Pandas), covering common data scenarios in bioinformatics",wordCount:314},{id:8,no:8,title:"R Language Primer: Data Structures, Vectorization, and Functions",date:"2026-08-04",category:["Programming","r"],tags:["R Language","Getting Started","Tutorial"],preview:"Core R fundamentals: environment setup, vectors/matrices/data frames/lists, vectorized operations, control flow, and function definitions",wordCount:256},{id:9,no:9,title:"R ggplot2 Data Visualization: Grammar of Graphics and Common Charts",date:"2026-08-04",category:["Programming","r"],tags:["R Language","ggplot2","Data Visualization","Tutorial"],preview:"Introduction to the ggplot2 Grammar of Graphics: scatter plots, boxplots, histograms, bar charts, and theme customization for publication-ready figures",wordCount:157},{id:10,no:10,title:"R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning",date:"2026-08-04",category:["Programming","r"],tags:["R Language","tidyverse","dplyr","Tutorial"],preview:"Chain dplyr verbs (filter/select/mutate/group_by/summarise/join) with pipes, and use tidyr to complete a full data cleaning workflow.",wordCount:245},{id:11,no:11,title:"Review of Cryo-EM Technology",date:"2026-08-04",category:["Structural Biology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Review"],preview:"The technological revolution of single-particle cryo-EM: direct electron detectors, the stability revolution, and the advent of the AI era, along with the future landscape of structural biology",wordCount:977},{id:12,no:12,title:"Cryo-EM Single Particle Analysis: Full Data Processing Workflow",date:"2026-08-04",category:["Structural Biology","cryoem"],tags:["Cryo-Electron Microscopy","cryo-EM","Data Processing","RELION"],preview:"From micrographs to atomic models: the complete workflow of cryo-EM single particle analysis (SPA), covering preprocessing, particle picking, 2D/3D classification, refinement, and resolution assessment",wordCount:765},{id:13,no:13,title:"Atomic Modeling and Refinement",date:"2026-08-04",category:["Structural Biology","modeling"],tags:["Atomic Modeling","Coot","phenix","Structure Refinement","Tutorial"],preview:"Complete workflow for atomic model building and refinement from cryo-EM density maps: automatic modeling, Coot manual correction, phenix refinement, model quality assessment metrics",wordCount:566},{id:14,no:14,title:"Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice",date:"2026-08-04",category:["Structural Biology","visualization"],tags:["Structure Visualization","PyMOL","ChimeraX","Tutorial"],preview:"Visualizing protein/nucleic acid structures with PyMOL and UCSF ChimeraX: PDB data retrieval, rendering modes, script batch processing, and cryo-EM density map display",wordCount:362}],Xb=Object.freeze(Object.defineProperty({__proto__:null,default:Kb},Symbol.toStringTag,{value:"Module"})),Yb=[{id:1,no:1,title:"虚拟筛选主流工具",date:"2026-08-04",category:["计算生物学","virtual-screening"],tags:["虚拟筛选","分子对接","AutoDock Vina","计算生物学"],preview:"虚拟筛选工具全景：AutoDock Vina、Glide、DOCK 等对接软件，受体/配体准备、打分函数、高通量筛选流程与验证策略",wordCount:897},{id:2,no:2,title:"Bash 编程：变量、条件、循环与实用脚本",date:"2026-08-04",category:["编程","bash"],tags:["Bash","Shell","脚本编程","教程"],preview:"Bash 脚本编程完整教程：变量与参数、条件判断、循环、函数、数组与字符串处理，附批处理与生信实用脚本",wordCount:259},{id:3,no:3,title:"Linux 命令行基础：文件系统、权限与文本处理",date:"2026-08-04",category:["编程","linux"],tags:["Linux","命令行","教程"],preview:"Linux 命令行核心技能：文件系统导航、文件操作、权限管理、文本处理三剑客（grep/sed/awk）与进程管理",wordCount:328},{id:4,no:4,title:"Python 进阶：函数、类与模块",date:"2026-08-04",category:["编程","python"],tags:["Python","进阶","教程"],preview:"深入 Python 函数式与面向对象编程：参数传递、lambda、装饰器、类与继承、异常处理、模块与包",wordCount:376},{id:5,no:5,title:"Python 编程入门：环境、语法与数据类型",date:"2026-08-04",category:["编程","python"],tags:["Python","入门","教程"],preview:"从零开始学习 Python：环境搭建、基础语法、内置数据类型与流程控制，附带可运行的示例代码",wordCount:495},{id:6,no:6,title:"Python 数据处理实战：文件 IO、正则与 NumPy/Pandas",date:"2026-08-04",category:["编程","python"],tags:["Python","NumPy","Pandas","数据处理","教程"],preview:"文件读写、正则表达式与科学计算三件套（NumPy/Pandas）的完整实战教程，覆盖生物信息学常见数据场景",wordCount:473},{id:7,no:7,title:"R 语言入门：数据结构、向量化与函数",date:"2026-08-04",category:["编程","r"],tags:["R语言","入门","教程"],preview:"R 语言核心基础：环境配置、向量/矩阵/数据框/列表、向量化运算、控制流与函数定义",wordCount:403},{id:8,no:8,title:"R ggplot2 数据可视化：图层语法与常用图表",date:"2026-08-04",category:["编程","r"],tags:["R语言","ggplot2","数据可视化","教程"],preview:"ggplot2 图层语法入门：散点图、箱线图、直方图、柱状图与主题美化，输出出版级图表",wordCount:245},{id:9,no:9,title:"R tidyverse 数据操作：dplyr 管道与数据清洗",date:"2026-08-04",category:["编程","r"],tags:["R语言","tidyverse","dplyr","教程"],preview:"用管道符串联 dplyr 动词（filter/select/mutate/group_by/summarise/join），配合 tidyr 完成完整的数据清洗流程",wordCount:324},{id:10,no:10,title:"冷冻电镜技术综述",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","综述"],preview:"冷冻电镜单颗粒分析的技术革命：直接电子探测器、稳定性革命与 AI 时代的到来，以及结构生物学的未来图景",wordCount:1485},{id:11,no:11,title:"冷冻电镜单颗粒分析：数据处理全流程",date:"2026-08-04",category:["结构生物学","cryoem"],tags:["冷冻电镜","cryo-EM","数据处理","RELION"],preview:"从显微照片到原子模型：冷冻电镜单颗粒分析（SPA）完整流程，涵盖预处理、颗粒挑选、2D/3D 分类、精修与分辨率评估",wordCount:1111},{id:12,no:12,title:"原子建模与精修",date:"2026-08-04",category:["结构生物学","modeling"],tags:["原子建模","Coot","phenix","结构精修","教程"],preview:"冷冻电镜密度图的原子模型搭建与精修全流程：自动建模、Coot 手动修正、phenix 精修、模型质量评估指标",wordCount:835},{id:13,no:13,title:"生物大分子结构可视化：PyMOL 与 ChimeraX 实战",date:"2026-08-04",category:["结构生物学","visualization"],tags:["结构可视化","PyMOL","ChimeraX","教程"],preview:"用 PyMOL 与 UCSF ChimeraX 可视化蛋白质/核酸结构：PDB 数据获取、渲染模式、脚本批处理与冷冻电镜密度图展示",wordCount:529},{id:14,no:14,title:"蛋白质设计主流工具",date:"",category:["计算生物学","protein-design"],tags:["蛋白质设计","Rosetta","AlphaFold","RFdiffusion","计算生物学"],preview:"蛋白质设计工具全景：Rosetta 能量优化、AlphaFold 结构预测、RFdiffusion 生成式设计、ESMFold 快速预测与 binder 设计工作流",wordCount:1002}],Qb=Object.freeze(Object.defineProperty({__proto__:null,default:Yb},Symbol.toStringTag,{value:"Module"})),Jb=[{name:"Plotit",title:"plotit",desc:"Declarative plotting package for R — a verb-prefix API based on ggplot2, in early development stage",date:"2025-12-30",tags:["R","ggplot2","Data Visualization","Declarative API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"Personal blog system — built with Vue 3 + Vite, supporting bilingual Chinese/English, theme switching, and Markdown article management",date:"2025-08-29",tags:["Vue","Vite","Blog","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],Zb=Object.freeze(Object.defineProperty({__proto__:null,default:Jb},Symbol.toStringTag,{value:"Module"})),s_=[{name:"Plotit",title:"plotit",desc:"R 声明式绘图包——基于 ggplot2 的动词前缀 API，早期开发阶段",date:"2025-12-30",tags:["R","ggplot2","数据可视化","声明式API"],github:"https://github.com/zorrooz/plotit",url:"https://zorrooz.github.io/plotit",status:"active",language:"R",stars:0,license:"MIT",version:""},{name:"ZorroozBlog",title:"zorrooz.github.io",desc:"个人博客系统——基于 Vue 3 + Vite 构建，支持中英双语、主题切换与 Markdown 文章管理",date:"2025-08-29",tags:["Vue","Vite","博客","Markdown"],github:"https://github.com/zorrooz/zorrooz.github.io",url:"https://zorrooz.github.io/",status:"active",language:"Vue",stars:0,license:"MIT",version:""}],n_=Object.freeze(Object.defineProperty({__proto__:null,default:s_},Symbol.toStringTag,{value:"Module"})),a_=JSON.parse(`[{"title":"Data Analysis","children":[{"title":"Data Exploration","children":[{"title":"R Language","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"A suite of R packages for data science, covering data import, cleaning, transformation, and visualization"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"High-performance data manipulation package, extremely fast for processing data frames with millions of rows"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"A tool for fast reading of CSV/TSV and other tabular files, with automatic column type inference"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"The fundamental library for scientific computing in Python, providing multidimensional arrays and vectorized operations"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"A powerful data analysis library based on Python, providing efficient data manipulation and processing capabilities"},{"name":"Polars","url":"https://www.pola.rs/","desc":"A high-performance data processing library based on Apache Arrow, supporting multiple language interfaces"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"An interactive notebook environment that supports mixing code, text, and visualizations"}]}]},{"title":"Data Visualization","children":[{"title":"R Language","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"A powerful data visualization package in R, based on the grammar of graphics"},{"name":"plotly (R)","url":"https://plotly.com/r/","desc":"The R interface to the interactive chart library, capable of generating web-interactive graphics"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"A multi-panel composition tool that easily combines ggplot graphics using + and / symbols"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"Publication-quality color schemes, providing carefully designed discrete palettes"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"A widely used plotting library in Python, supporting many chart types"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"A statistical visualization library based on Matplotlib, with built-in attractive themes and statistical plots"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"An interactive visualization library supporting scatter, 3D, maps, and many other chart types"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"A declarative statistical visualization library based on Vega-Lite"}]}]},{"title":"Statistical Analysis","children":[{"title":"R Language","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"The built-in statistical analysis functionality in R, covering a wide range of statistical methods and models"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"Pipe-friendly wrappers for statistical tests (t-tests, ANOVA, rank-sum tests, etc.)"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"A tool for tidying statistical model outputs into clean data frames"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"Linear and nonlinear mixed-effects models, suitable for repeated-measures data"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"The statistics module in the Python scientific computing library SciPy, providing many statistical distributions and test methods"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"A statistical modeling library supporting regression, time series, hypothesis testing, and more"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"A simple and user-friendly statistical testing library covering common parametric and nonparametric tests"}]}]},{"title":"Machine Learning","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"A machine learning library based on Python, providing a rich set of algorithms and tools"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"A deep learning framework supporting dynamic computation graphs and efficient tensor operations"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"A deep learning framework with a mature ecosystem, supporting production-level deployment"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"A gradient boosting tree library, the top choice for tabular data competitions and engineering"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"A gradient boosting framework from Microsoft, with fast training speed and low memory usage"}]},{"title":"R Language","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"A unified modeling framework in R, covering data preprocessing, modeling, and evaluation"}]}]}]},{"title":"Omics","children":[{"title":"Data Platforms","children":[{"title":"Sequence Databases","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"The National Center for Biotechnology Information in the United States, providing important databases such as GenBank"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"The European Molecular Biology Laboratory, providing a variety of bioinformatics resources and tools"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"The DNA Data Bank of Japan, providing storage and access for nucleic acid and protein sequence data"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"The most comprehensive database for protein sequences and functional annotation"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"A genome annotation and comparative genomics database for vertebrates"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"A genome browser supporting visualization of multiple species genomes and custom tracks"}]},{"title":"Sequencing Data","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"The NCBI Sequence Read Archive, storing raw high-throughput sequencing data"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"A gene expression database containing microarray and sequencing expression data"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"The European Nucleotide Archive, Europe's storage center for raw sequencing data"}]},{"title":"Protein Interactions and Pathways","items":[{"name":"STRING","url":"https://string-db.org/","desc":"A protein-protein interaction database supporting multiple species"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"A comprehensive database of pathways, genes, and compounds"},{"name":"GO","url":"http://geneontology.org/","desc":"Gene Ontology, providing a standardized system for gene function annotation"}]}]},{"title":"Analysis Software","children":[{"title":"Command Line","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"A quality assessment tool for sequencing data, generating visual QC reports"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"A read trimming tool for sequencing data, removing adapters and low-quality bases"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"A command-line adapter removal tool supporting multiple sequencing platforms"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"A short-read alignment tool, the standard choice for genome alignment"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"An ultra-fast splice-aware aligner for RNA-seq"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"A graph-based index alignment tool for RNA-seq"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"A toolset for processing SAM/BAM files, the standard for post-alignment processing"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"A tool for processing VCF/BCF variant files, supporting filtering and statistics"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"The Genome Analysis Toolkit, the industry standard for variant detection"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"Aggregates multiple QC reports into a single interactive HTML report"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"The core Python library for bioinformatics, covering sequences, structures, and file parsing"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"A Python library for bioinformatics, providing various tools for biological data processing and analysis"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"A Python library for single-cell transcriptomics analysis, deeply integrated with AnnData"}]},{"title":"R Language","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"A collection of bioinformatics software packages based on R, covering genomics, transcriptomics, and other fields"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"A mainstream R package for single-cell transcriptomics analysis, providing a complete analysis workflow"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"The standard tool for differential expression analysis based on the negative binomial model"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"A differential expression analysis tool based on empirical Bayes"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"A differential analysis tool using linear models, applicable to both microarray and sequencing data"}]}]}]},{"title":"Cryo-EM Structure Determination","children":[{"title":"Software Tools","children":[{"title":"3D Reconstruction","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"A Bayesian approach-based software for cryo-EM 3D reconstruction, widely used for high-resolution structure determination"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"A high-performance cryo-EM data processing platform with an intuitive user interface"},{"name":"cisTEM","url":"https://cistem.org/","desc":"A free, open-source cryo-EM processing software supporting a complete single-particle workflow"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"A modular cryo-EM processing suite supporting high-resolution reconstruction"}]},{"title":"Preprocessing","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"A frame-by-frame motion correction tool that eliminates sample drift caused by the electron beam"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"A classic tool for CTF estimation, evaluating defocus and astigmatism"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"A GPU-accelerated CTF estimation tool with high speed and precision"},{"name":"Warp","url":"http://www.warpem.com/","desc":"A fast preprocessing tool integrating motion correction, CTF estimation, and particle picking"}]},{"title":"Particle Picking","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"A deep learning-based particle picking tool supporting real-time picking"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"A CNN-based particle picking tool with excellent performance on low signal-to-noise ratio samples"}]},{"title":"Modeling and Visualization","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"A molecular visualization software, the standard tool for scripted rendering and structure analysis"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"A modern molecular visualization software with outstanding density map processing capabilities"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"A manual model building tool for interactive adjustment of density maps and models"},{"name":"phenix","url":"https://phenix-online.org/","desc":"A structure refinement and validation suite supporting cryo-EM and crystallography"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"An interactive real-time molecular dynamics modeling tool based on ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"A deep learning-based automatic atomic modeling tool for density maps"}]}]},{"title":"Database Resources","children":[{"title":"Structures and Density Maps","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"The Electron Microscopy Data Bank, storing and distributing 3D density maps from cryo-EM"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"The Electron Microscopy Public Image Archive, providing raw movies and particle data"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"The Protein Data Bank, containing biological macromolecular structures determined by X-ray, NMR, and cryo-EM"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"The worldwide PDB coordination body, uniformly managing structural data standards"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"The European PDB node, providing structural data and visualization tools"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"The AlphaFold protein structure database, covering hundreds of millions of proteins"}]}]}]}]`),e_=Object.freeze(Object.defineProperty({__proto__:null,default:a_},Symbol.toStringTag,{value:"Module"})),t_=JSON.parse('[{"title":"数据分析","children":[{"title":"数据探索","children":[{"title":"R 语言","items":[{"name":"tidyverse","url":"https://www.tidyverse.org/","desc":"R 语言中用于数据科学的一整套包，涵盖数据导入、清洗、转换和可视化"},{"name":"data.table","url":"https://rdatatable.gitlab.io/data.table/","desc":"高性能数据操作包，百万行级数据框处理速度极快"},{"name":"readr","url":"https://readr.tidyverse.org/","desc":"快速读取 CSV/TSV 等表格文件的工具，自动推断列类型"}]},{"title":"Python","items":[{"name":"NumPy","url":"https://numpy.org/","desc":"Python 科学计算基础库，提供多维数组与向量化运算"},{"name":"Pandas","url":"https://pandas.pydata.org/","desc":"基于 Python 的强大数据分析库，提供高效的数据操作和处理功能"},{"name":"Polars","url":"https://www.pola.rs/","desc":"基于 Apache Arrow 的高性能数据处理库，支持多语言接口"},{"name":"Jupyter","url":"https://jupyter.org/","desc":"交互式笔记本环境，支持代码、文本与可视化混排"}]}]},{"title":"数据可视化","children":[{"title":"R 语言","items":[{"name":"ggplot2","url":"https://ggplot2.tidyverse.org/","desc":"R 语言中功能强大的数据可视化包，基于语法图形学"},{"name":"plotly（R）","url":"https://plotly.com/r/","desc":"交互式图表库的 R 接口，可生成网页端可交互图形"},{"name":"patchwork","url":"https://patchwork.data-imaginist.com/","desc":"多图拼接工具，用 + 和 / 符号轻松组合 ggplot 图形"},{"name":"RColorBrewer","url":"https://cran.r-project.org/web/packages/RColorBrewer/","desc":"出版级配色方案，提供精心设计的离散色板"}]},{"title":"Python","items":[{"name":"Matplotlib","url":"https://matplotlib.org/","desc":"Python 中广泛使用的绘图库，支持多种图表类型"},{"name":"Seaborn","url":"https://seaborn.pydata.org/","desc":"基于 Matplotlib 的统计可视化库，内置美观主题与统计图"},{"name":"Plotly","url":"https://plotly.com/python/","desc":"交互式可视化库，支持散点、3D、地图等多种图表"},{"name":"Altair","url":"https://altair-viz.github.io/","desc":"基于 Vega-Lite 的声明式统计可视化库"}]}]},{"title":"统计分析","children":[{"title":"R 语言","items":[{"name":"R Stats","url":"https://www.r-project.org/","desc":"R 语言内置的统计分析功能，涵盖广泛的统计方法和模型"},{"name":"rstatix","url":"https://rpkgs.datanovia.com/rstatix/","desc":"管道友好的统计检验封装（t 检验、方差分析、秩和检验等）"},{"name":"broom","url":"https://broom.tidymodels.org/","desc":"把统计模型输出整理为整洁数据框的工具"},{"name":"lme4","url":"https://cran.r-project.org/web/packages/lme4/","desc":"线性与非线性混合效应模型，适用于重复测量数据"}]},{"title":"Python","items":[{"name":"SciPy Stats","url":"https://scipy.org/","desc":"Python 科学计算库 SciPy 中的统计模块，提供多种统计分布和检验方法"},{"name":"statsmodels","url":"https://www.statsmodels.org/","desc":"统计建模库，支持回归、时间序列、假设检验等"},{"name":"pingouin","url":"https://pingouin-stats.org/","desc":"简洁易用的统计检验库，覆盖常用参数与非参数检验"}]}]},{"title":"机器学习","children":[{"title":"Python","items":[{"name":"scikit-learn","url":"https://scikit-learn.org/","desc":"基于 Python 的机器学习库，提供丰富的算法和工具"},{"name":"PyTorch","url":"https://pytorch.org/","desc":"深度学习框架，支持动态计算图和高效的张量运算"},{"name":"TensorFlow","url":"https://www.tensorflow.org/","desc":"深度学习框架，生态成熟，支持生产级部署"},{"name":"XGBoost","url":"https://xgboost.readthedocs.io/","desc":"梯度提升树库，表格数据竞赛与工程的首选"},{"name":"LightGBM","url":"https://lightgbm.readthedocs.io/","desc":"微软出品的梯度提升框架，训练速度快、内存占用低"}]},{"title":"R 语言","items":[{"name":"tidymodels","url":"https://www.tidymodels.org/","desc":"R 语言统一的建模框架，覆盖数据预处理、建模与评估"}]}]}]},{"title":"组学","children":[{"title":"数据平台","children":[{"title":"序列数据库","items":[{"name":"NCBI","url":"https://www.ncbi.nlm.nih.gov/","desc":"美国国家生物技术信息中心，提供 GenBank 等重要数据库"},{"name":"EMBL-EBI","url":"https://www.ebi.ac.uk/","desc":"欧洲分子生物学实验室，提供多种生物信息学资源和工具"},{"name":"DDBJ","url":"https://www.ddbj.nig.ac.jp/","desc":"日本 DNA 数据库，提供核酸、蛋白序列数据的存储和访问"},{"name":"UniProt","url":"https://www.uniprot.org/","desc":"最全面的蛋白质序列与功能注释数据库"},{"name":"Ensembl","url":"https://www.ensembl.org/","desc":"脊椎动物基因组注释与比较基因组学数据库"},{"name":"UCSC Genome Browser","url":"https://genome.ucsc.edu/","desc":"基因组浏览器，支持多物种基因组可视化与自定义轨道"}]},{"title":"测序数据","items":[{"name":"SRA","url":"https://www.ncbi.nlm.nih.gov/sra","desc":"NCBI 原始测序数据仓库，存储高通量测序原始数据"},{"name":"GEO","url":"https://www.ncbi.nlm.nih.gov/geo/","desc":"基因表达数据库，收录芯片与测序表达数据"},{"name":"ENA","url":"https://www.ebi.ac.uk/ena","desc":"欧洲核酸档案库，欧洲的原始测序数据存储中心"}]},{"title":"蛋白互作与通路","items":[{"name":"STRING","url":"https://string-db.org/","desc":"蛋白-蛋白相互作用数据库，支持多种物种"},{"name":"KEGG","url":"https://www.kegg.jp/","desc":"通路、基因与化合物综合数据库"},{"name":"GO","url":"http://geneontology.org/","desc":"基因本体论，提供标准化的基因功能注释体系"}]}]},{"title":"分析软件","children":[{"title":"命令行","items":[{"name":"FastQC","url":"https://www.bioinformatics.babraham.ac.uk/projects/fastqc/","desc":"测序数据质量评估工具，生成可视化质控报告"},{"name":"Trimmomatic","url":"http://www.usadellab.org/cms/?page=trimmomatic","desc":"测序读段修剪工具，去除接头与低质量碱基"},{"name":"cutadapt","url":"https://cutadapt.readthedocs.io/","desc":"命令行接头去除工具，支持多种测序平台"},{"name":"BWA","url":"https://github.com/lh3/bwa","desc":"短读段比对工具，基因组比对的标准选择"},{"name":"STAR","url":"https://github.com/alexdobin/STAR","desc":"RNA-seq 剪接感知比对器，速度极快"},{"name":"HISAT2","url":"https://daehwankimlab.github.io/hisat2/","desc":"基于图索引的 RNA-seq 比对工具"},{"name":"SAMtools","url":"https://www.htslib.org/","desc":"SAM/BAM 文件处理工具集，比对后处理的标准"},{"name":"bcftools","url":"https://samtools.github.io/bcftools/","desc":"VCF/BCF 变异文件处理工具，支持过滤与统计"},{"name":"GATK","url":"https://gatk.broadinstitute.org/","desc":"基因组分析工具包，变异检测的行业标准"},{"name":"MultiQC","url":"https://multiqc.info/","desc":"把多个质控报告汇总为一个交互式 HTML 报告"}]},{"title":"Python","items":[{"name":"Biopython","url":"https://biopython.org/","desc":"生物信息学核心 Python 库，涵盖序列、结构与文件解析"},{"name":"scikit-bio","url":"https://scikit-bio.org/","desc":"用于生物信息学的 Python 库，提供多种生物数据处理和分析工具"},{"name":"scanpy","url":"https://scanpy.readthedocs.io/","desc":"单细胞转录组分析 Python 库，与 AnnData 深度集成"}]},{"title":"R 语言","items":[{"name":"Bioconductor","url":"https://bioconductor.org/","desc":"基于 R 语言的生物信息学软件包集合，涵盖基因组学、转录组学等领域"},{"name":"Seurat","url":"https://satijalab.org/seurat/","desc":"单细胞转录组分析主流 R 包，提供完整分析流程"},{"name":"DESeq2","url":"https://bioconductor.org/packages/release/bioc/html/DESeq2.html","desc":"基于负二项模型的差异表达分析标准工具"},{"name":"edgeR","url":"https://bioconductor.org/packages/release/bioc/html/edgeR.html","desc":"基于经验贝叶斯的差异表达分析工具"},{"name":"limma","url":"https://bioconductor.org/packages/release/bioc/html/limma.html","desc":"线性模型差异分析工具，芯片与测序数据通用"}]}]}]},{"title":"电镜结构解析","children":[{"title":"软件工具","children":[{"title":"三维重构","items":[{"name":"RELION","url":"https://www3.mrc-lmb.cam.ac.uk/groups/scheres/index.php","desc":"基于贝叶斯方法的冷冻电镜三维重构软件，广泛应用于高分辨率结构解析"},{"name":"cryoSPARC","url":"https://cryoem.slac.stanford.edu/","desc":"高性能的冷冻电镜数据处理平台，提供直观的用户界面"},{"name":"cisTEM","url":"https://cistem.org/","desc":"免费开源的冷冻电镜处理软件，支持完整的单颗粒流程"},{"name":"SPHIRE","url":"https://sphire.mpg.de/","desc":"模块化冷冻电镜处理套件，支持高分辨率重构"}]},{"title":"预处理","items":[{"name":"MotionCor2","url":"https://emcore.ucsf.edu/ucsf-software","desc":"逐帧运动校正工具，消除电子束导致的样品漂移"},{"name":"CTFFIND4","url":"https://grigoriefflab.umassmed.edu/ctffind4","desc":"CTF 估计经典工具，评估欠焦与像散"},{"name":"Gctf","url":"https://www.mrc-lmb.cam.ac.uk/kzhang/Gctf/","desc":"GPU 加速的 CTF 估计工具，速度快精度高"},{"name":"Warp","url":"http://www.warpem.com/","desc":"集运动校正、CTF 与颗粒挑选于一体的快速预处理工具"}]},{"title":"颗粒挑选","items":[{"name":"crYOLO","url":"https://cryolo.readthedocs.io/","desc":"基于深度学习的颗粒挑选工具，支持实时挑选"},{"name":"Topaz","url":"https://topaz.readthedocs.io/","desc":"基于 CNN 的颗粒挑选工具，对低信噪比样品表现优异"}]},{"title":"建模与可视化","items":[{"name":"PyMOL","url":"https://pymol.org/","desc":"分子可视化软件，脚本化渲染与结构分析的标准工具"},{"name":"UCSF ChimeraX","url":"https://www.cgl.ucsf.edu/chimerax/","desc":"现代分子可视化软件，密度图处理能力突出"},{"name":"Coot","url":"https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/","desc":"手动模型修正工具，密度图与模型交互调整"},{"name":"phenix","url":"https://phenix-online.org/","desc":"结构精修与验证套件，支持冷冻电镜与晶体学"},{"name":"ISOLDE","url":"https://isolde.cimr.cam.ac.uk/","desc":"交互式实时分子动力学建模工具，基于 ChimeraX"},{"name":"ModelAngelo","url":"https://github.com/3dem/model-angelo","desc":"基于深度学习的密度图自动原子建模工具"}]}]},{"title":"数据库资源","children":[{"title":"结构与密度图","items":[{"name":"EMDB","url":"https://www.ebi.ac.uk/emdb/","desc":"电子显微镜数据银行，存储和分发冷冻电镜三维密度图"},{"name":"EMPIAR","url":"https://www.ebi.ac.uk/empiar/","desc":"冷冻电镜原始数据档案，提供原始电影与颗粒数据"},{"name":"RCSB PDB","url":"https://www.rcsb.org/","desc":"蛋白质数据银行，包含通过X射线、NMR和冷冻电镜确定的生物大分子结构"},{"name":"wwPDB","url":"https://www.wwpdb.org/","desc":"全球 PDB 协调机构，统一管理结构数据标准"},{"name":"PDBe","url":"https://www.ebi.ac.uk/pdbe/","desc":"欧洲 PDB 节点，提供结构数据与可视化工具"},{"name":"AlphaFold DB","url":"https://alphafold.ebi.ac.uk/","desc":"AlphaFold 预测蛋白质结构数据库，覆盖数亿蛋白质"}]}]}]}]'),l_=Object.freeze(Object.defineProperty({__proto__:null,default:t_},Symbol.toStringTag,{value:"Module"})),o_=[{name:"Advanced",count:1},{name:"AlphaFold",count:1},{name:"Atomic Modeling",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Command Line",count:1},{name:"Computational Biology",count:2},{name:"Coot",count:1},{name:"Cryo-Electron Microscopy",count:2},{name:"cryo-EM",count:2},{name:"Data Processing",count:2},{name:"Data Visualization",count:1},{name:"dplyr",count:1},{name:"Getting Started",count:2},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"Molecular Docking",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"Protein Design",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R Language",count:3},{name:"RELION",count:1},{name:"Review",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Scripting",count:1},{name:"Shell",count:1},{name:"Structure Refinement",count:1},{name:"Structure Visualization",count:1},{name:"tidyverse",count:1},{name:"Tutorial",count:10},{name:"Virtual Screening",count:1}],p_=Object.freeze(Object.defineProperty({__proto__:null,default:o_},Symbol.toStringTag,{value:"Module"})),c_=[{name:"蛋白质设计",count:1},{name:"分子对接",count:1},{name:"计算生物学",count:2},{name:"脚本编程",count:1},{name:"教程",count:10},{name:"结构精修",count:1},{name:"结构可视化",count:1},{name:"进阶",count:1},{name:"冷冻电镜",count:2},{name:"命令行",count:1},{name:"入门",count:2},{name:"数据处理",count:2},{name:"数据可视化",count:1},{name:"虚拟筛选",count:1},{name:"原子建模",count:1},{name:"综述",count:1},{name:"AlphaFold",count:1},{name:"AutoDock Vina",count:1},{name:"Bash",count:1},{name:"ChimeraX",count:1},{name:"Coot",count:1},{name:"cryo-EM",count:2},{name:"dplyr",count:1},{name:"ggplot2",count:1},{name:"Linux",count:1},{name:"NumPy",count:1},{name:"Pandas",count:1},{name:"phenix",count:1},{name:"PyMOL",count:1},{name:"Python",count:3},{name:"R语言",count:3},{name:"RELION",count:1},{name:"RFdiffusion",count:1},{name:"Rosetta",count:1},{name:"Shell",count:1},{name:"tidyverse",count:1}],r_=Object.freeze(Object.defineProperty({__proto__:null,default:c_},Symbol.toStringTag,{value:"Module"})),i_=[{name:"StructureOfProteinDemo",title:"Protein Structure Determination Demo Project",desc:"Protein structure placeholder demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],u_=Object.freeze(Object.defineProperty({__proto__:null,default:i_},Symbol.toStringTag,{value:"Module"})),h_=[{name:"StructureOfProteinDemo",title:"蛋白质结构解析演示课题",desc:"蛋白质结构占位demo",date:"2026-04-12",tags:["Protein","Structure","Cryo-EM"],doi:"10.1234/demo-structure.2026",url:"https://www.demo-structure.org",status:"completed",journal:"Journal of Demo Science",year:2026,authors:["Demo Author A","Demo Author B"]}],d_=Object.freeze(Object.defineProperty({__proto__:null,default:h_},Symbol.toStringTag,{value:"Module"})),m_=`<h1>Mainstream Tools for Protein Design</h1>
<p>Protein design is an important direction in computational biology: given a target function, design an amino acid sequence that folds into a specific structure. This article reviews the mainstream tool spectrum from classical energy functions to generative AI.</p>
<h2>1. Two Paradigms of Design Problems</h2>
<h3>1.1 Sequence Design</h3>
<p>Given a target structure (backbone), find the most stable amino acid sequence:</p>
<pre><code>structure (backbone) → sequence design → amino acid sequence → experimental validation
</code></pre>
<h3>1.2 Structure Design (De Novo Design)</h3>
<p>Design new structures/functions from scratch:</p>
<pre><code>functional requirement → structure generation → sequence design → experimental validation
</code></pre>
<h2>2. Rosetta: The Energy Function Paradigm</h2>
<p>Rosetta is the "veteran empire" in structure prediction and design, based on a hybrid energy function of physics and knowledge.</p>
<h3>2.1 Core Energy Function</h3>
<ul>
<li><strong>REF2015</strong>: the classic all-atom energy function (van der Waals, hydrogen bonding, solvation, electrostatics, and dihedral angle preferences)</li>
<li>Score = weighted sum of energy terms; lower is better</li>
<li>Rosetta's "design" is essentially <strong>energy minimization</strong>: searching for low-energy states in sequence/structure space</li>
</ul>
<h3>2.2 Common Protocols</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Sequence design (FixBB: fixed backbone design)</span>
rosetta_scripts.default.linuxgccrelease \\
  -s input.pdb \\
  -parser:protocol fixbb.xml \\
  -out:prefix designed_

<span class="hljs-comment"># Core snippet of fixbb.xml</span>
<span class="hljs-comment"># &#x3C;TaskOperations></span>
<span class="hljs-comment">#   &#x3C;RestrictToRepacking name="repack"/></span>
<span class="hljs-comment">#   &#x3C;OperateOnResidueSubset name="no_design_core" selector="..." /></span>
<span class="hljs-comment"># &#x3C;/TaskOperations></span>
<span class="hljs-comment"># &#x3C;RosettaScripts></span>
<span class="hljs-comment">#   &#x3C;PackRotamersMover name="pack" task_operations="repack"/></span>
<span class="hljs-comment"># &#x3C;/RosettaScripts></span>
</code></pre>
<h3>2.3 Classic Applications of Rosetta</h3>
<ul>
<li><strong>Enzyme design</strong>: design of scaffolds for catalytic active sites (e.g., Kemp eliminase)</li>
<li><strong>Protein-protein interface design</strong>: engineering binding specificity</li>
<li><strong>Binder design</strong>: high-affinity binding proteins</li>
<li><strong>Thermostability engineering</strong>: computational mutations to increase Tm</li>
</ul>
<h3>2.4 Advantages and Limitations</h3>
<table>
<thead>
<tr>
<th>Advantages</th>
<th>Limitations</th>
</tr>
</thead>
<tbody>
<tr>
<td>Physically interpretable, finely controllable</td>
<td>Computationally intensive, requires empirical parameter tuning</td>
</tr>
<tr>
<td>Supports arbitrary backbones and modifications</td>
<td>Limited prediction for flexible regions</td>
</tr>
<tr>
<td>Mature ecosystem (ample documentation, tutorials)</td>
<td>Steep learning curve</td>
</tr>
</tbody>
</table>
<h2>3. AlphaFold: Prediction as a Design Aid</h2>
<h3>3.1 AlphaFold2 (2021)</h3>
<ul>
<li>End-to-end deep learning: sequence → structure (MSA + attention mechanism)</li>
<li>Accuracy: atomic-level accuracy at CASP14 (GDT > 90)</li>
<li>Value for design: <strong>rapidly verify the foldability of designed sequences</strong></li>
</ul>
<h3>3.2 AlphaFold3 (2024)</h3>
<ul>
<li>Unified framework for predicting protein, nucleic acid, small molecule, and ion complexes</li>
<li>Introduces diffusion models and pairwise representations</li>
<li>Significance: designed targets (binder–ligand complex structures) can be directly predicted</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ColabFold: a lightweight implementation of AlphaFold2 (GPU-friendly)</span>
colabfold_batch input.fasta out_dir \\
  --num-recycle 3 --model-type alphafold2_multimer_v3
</code></pre>
<h3>3.3 Key Usage: Validating Designs</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Confidence indicators for designed sequences</span>
<span class="hljs-comment"># pLDDT > 90: high confidence</span>
<span class="hljs-comment"># PAE &#x3C; 5: inter-domain relative positions trustworthy</span>
<span class="hljs-comment"># Design validation pipeline: sequence → AF2 prediction → comparison with target structure (RMSD)</span>
</code></pre>
<h2>4. Generative Design: RFdiffusion and ProteinMPNN</h2>
<h3>4.1 RFdiffusion (2023, Baker Lab)</h3>
<p>Structure generator based on <strong>denoising diffusion models</strong>:</p>
<ul>
<li>Input: functional motifs, shape constraints, binding targets</li>
<li>Output: novel protein backbones</li>
<li>Applications: binder design, symmetric oligomers, enzyme active-site scaffolds</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RFdiffusion example: designing a binding protein</span>
run_inference.py \\
  scaffoldguided.scaffoldguided=True \\
  scaff_loader.contigs_map.json=contigs.json \\
  inference.output_prefix=outputs/design \\
  inference.num_designs=100 \\
  potentials.guidance_scale=3.0
</code></pre>
<h3>4.2 ProteinMPNN (2022)</h3>
<p><strong>Sequence design neural network</strong>: given a backbone → sequence, forming a "generate-encode" loop with RFdiffusion:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ProteinMPNN sequence design</span>
python protein_mpnn_run.py \\
  --pdbpath scaffold.pdb \\
  --out_folder mpnn_out \\
  --num_seq_per_target 8 \\
  --batch_size 4
</code></pre>
<p><strong>Standard workflow (Baker Lab paradigm)</strong>:</p>
<pre><code>① RFdiffusion: generate a backbone
   ↓
② ProteinMPNN: backbone → multiple candidate sequences (sampling temperature 0.1-0.3)
   ↓
③ AF2 reverse validation: sequence → predicted structure → RMSD comparison with backbone
   ↓
④ Screen candidates with high pLDDT / low RMSD
   ↓
⑤ Experimental expression validation (yeast/phage display, ELISA)
</code></pre>
<h3>4.3 Other Generative Tools</h3>
<ul>
<li><strong>Chroma</strong>: another diffusion model (generation + sequence design in one)</li>
<li><strong>Genie</strong>: a faster protein diffusion model</li>
<li><strong>ESMFold</strong> (Meta, 2022): ultra-fast prediction without MSA (&#x3C;1 second per sequence), suitable for large-scale design screening</li>
<li><strong>ESM3</strong>: multimodal generation (sequence + structure + function)</li>
</ul>
<h2>5. Special Topic: Binder Design</h2>
<p>Binder (binding protein) design is currently the most active application scenario:</p>
<h3>5.1 Classic Workflow (Combining Existing Tools)</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Hotspot analysis on target surface (Rosetta)</span>
<span class="hljs-comment"># 2. Hotspot-guided binder backbone generation (RFdiffusion / or hotspot constraints)</span>
<span class="hljs-comment"># 3. ProteinMPNN sequence design</span>
<span class="hljs-comment"># 4. AF2 complex prediction validation (binder + target docking)</span>
<span class="hljs-comment"># 5. Experimental screening (yeast display, FACS sorting)</span>
</code></pre>
<h3>5.2 Key Considerations</h3>
<ul>
<li>Determination of binding hotspot residues: experimental mutation scanning (alanine scan) or computational (Rosetta ddG)</li>
<li>Binding affinity prediction: <code>FoldX</code>, <code>Rosetta interface_ddg</code>, AF2 PAE</li>
<li>Multiple "design-validate" iterations: experimental feedback returns to computation</li>
</ul>
<h2>6. Quick Reference Table for Tool Selection</h2>
<table>
<thead>
<tr>
<th>Task</th>
<th>Primary Tool</th>
<th>Alternative</th>
</tr>
</thead>
<tbody>
<tr>
<td>Structure prediction</td>
<td>AlphaFold2/3 (ColabFold)</td>
<td>ESMFold (speed)</td>
</tr>
<tr>
<td>Backbone generation</td>
<td>RFdiffusion</td>
<td>Chroma</td>
</tr>
<tr>
<td>Sequence design</td>
<td>ProteinMPNN</td>
<td>Rosetta FixBB</td>
</tr>
<tr>
<td>Thermostability design</td>
<td>Rosetta (ddg)</td>
<td>FoldX</td>
</tr>
<tr>
<td>Binding affinity assessment</td>
<td>Rosetta interface_ddg</td>
<td>FoldX / AF2 PAE</td>
</tr>
<tr>
<td>Interface design</td>
<td>Rosetta protocols</td>
<td>RFdiffusion + MPNN</td>
</tr>
<tr>
<td>Large-scale screening</td>
<td>ESMFold + ProteinMPNN</td>
<td>—</td>
</tr>
</tbody>
</table>
<h2>7. Experimental Validation Loop</h2>
<p>Computational design must return to experiments:</p>
<pre><code>Expression (E. coli / yeast) → Purification → Characterization
  → SEC / DSC (stability)
  → Surface plasmon resonance SPR / BLI (affinity)
  → Structural validation (crystallography / cryo-EM / AlphaFold comparison)
</code></pre>
<p><strong>Design success rate reference</strong>: Literature reports that AF2-guided RFdiffusion binder design can achieve experimental positive rates of 10–20% (far exceeding traditional methods), but every successful case involves multiple rounds of iteration.</p>
<h2>8. Summary</h2>
<ul>
<li><strong>Classic paradigm</strong>: Rosetta energy function (interpretable, adjustable)</li>
<li><strong>Prediction paradigm</strong>: AlphaFold family (validation, complex prediction)</li>
<li><strong>Generative paradigm</strong>: RFdiffusion + ProteinMPNN + AF2 validation loop (currently mainstream)</li>
<li>The ultimate measure of design success is always <strong>experimental validation</strong></li>
</ul>
<p>The next article will introduce mainstream tools for virtual screening: a complete toolchain for docking small molecules to targets.</p>`,j_=`<h1>蛋白质设计主流工具</h1>
<p>蛋白质设计（Protein Design）是计算生物学的重要方向：给定目标功能，设计出能折叠成特定结构的氨基酸序列。本文梳理从经典能量函数到生成式 AI 的主流工具谱系。</p>
<h2>1. 设计问题的两种范式</h2>
<h3>1.1 序列设计（Sequence Design）</h3>
<p>给定目标结构（骨架），寻找最稳定的氨基酸序列：</p>
<pre><code>结构（backbone）→ 序列设计 → 氨基酸序列 → 实验验证
</code></pre>
<h3>1.2 结构设计（De novo Design）</h3>
<p>从零设计新结构/新功能：</p>
<pre><code>功能需求 → 结构生成 → 序列设计 → 实验验证
</code></pre>
<h2>2. Rosetta：能量函数范式</h2>
<p>Rosetta 是结构预测与设计领域的"老牌帝国"，基于物理与知识混合的能量函数。</p>
<h3>2.1 核心能量函数</h3>
<ul>
<li><strong>REF2015</strong>：经典全原子能量函数（范德华、氢键、溶剂化、静电、二面角偏好）</li>
<li>打分 = 各项能量加权和，越低越好</li>
<li>Rosetta 的"设计"本质是<strong>能量最小化</strong>：在序列/结构空间搜索低能态</li>
</ul>
<h3>2.2 常用协议</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 序列设计（FixBB：固定骨架设计）</span>
rosetta_scripts.default.linuxgccrelease \\
  -s input.pdb \\
  -parser:protocol fixbb.xml \\
  -out:prefix designed_

<span class="hljs-comment"># fixbb.xml 核心片段</span>
<span class="hljs-comment"># &#x3C;TaskOperations></span>
<span class="hljs-comment">#   &#x3C;RestrictToRepacking name="repack"/></span>
<span class="hljs-comment">#   &#x3C;OperateOnResidueSubset name="no_design_core" selector="..." /></span>
<span class="hljs-comment"># &#x3C;/TaskOperations></span>
<span class="hljs-comment"># &#x3C;RosettaScripts></span>
<span class="hljs-comment">#   &#x3C;PackRotamersMover name="pack" task_operations="repack"/></span>
<span class="hljs-comment"># &#x3C;/RosettaScripts></span>
</code></pre>
<h3>2.3 Rosetta 的经典应用</h3>
<ul>
<li><strong>酶设计</strong>：催化活性位点支架设计（如 Kemp 消除酶）</li>
<li><strong>蛋白-蛋白界面设计</strong>：结合特异性改造</li>
<li><strong>binder 设计</strong>：高亲和力结合蛋白</li>
<li><strong>热稳定性改造</strong>：计算突变提高 Tm</li>
</ul>
<h3>2.4 优点与局限</h3>
<table>
<thead>
<tr>
<th>优点</th>
<th>局限</th>
</tr>
</thead>
<tbody>
<tr>
<td>物理可解释、可精细调控</td>
<td>计算量大、需要经验调参</td>
</tr>
<tr>
<td>支持任意骨架与修饰</td>
<td>对柔性区域预测有限</td>
</tr>
<tr>
<td>生态成熟（文档、教程多）</td>
<td>学习曲线陡峭</td>
</tr>
</tbody>
</table>
<h2>3. AlphaFold：预测即设计辅助</h2>
<h3>3.1 AlphaFold2（2021）</h3>
<ul>
<li>端到端深度学习：序列 → 结构（MSA + 注意力机制）</li>
<li>精度：CASP14 原子级精度（GDT > 90）</li>
<li>对设计的价值：<strong>快速验证设计序列的可折叠性</strong></li>
</ul>
<h3>3.2 AlphaFold3（2024）</h3>
<ul>
<li>统一框架预测蛋白质、核酸、小分子、离子复合物</li>
<li>引入扩散模型与配对表示</li>
<li>意义：设计靶标（binder 与配体的复合物结构）可直接预测</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ColabFold：AlphaFold2 的轻量实现（GPU 友好）</span>
colabfold_batch input.fasta out_dir \\
  --num-recycle 3 --model-type alphafold2_multimer_v3
</code></pre>
<h3>3.3 关键用法：验证设计</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 设计序列的可信度指标</span>
<span class="hljs-comment"># pLDDT > 90：高置信</span>
<span class="hljs-comment"># PAE &#x3C; 5：结构域间相对位置可信</span>
<span class="hljs-comment"># 设计验证流程：序列 → AF2 预测 → 与目标结构比对（RMSD）</span>
</code></pre>
<h2>4. 生成式设计：RFdiffusion 与 ProteinMPNN</h2>
<h3>4.1 RFdiffusion（2023，Baker Lab）</h3>
<p>基于<strong>去噪扩散模型</strong>的结构生成器：</p>
<ul>
<li>输入：功能基序（motif）、形状约束、结合靶点</li>
<li>输出：全新蛋白质骨架（backbone）</li>
<li>应用：binder 设计、对称寡聚体、酶活性位点支架</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RFdiffusion 示例：设计结合蛋白</span>
run_inference.py \\
  scaffoldguided.scaffoldguided=True \\
  scaff_loader.contigs_map.json=contigs.json \\
  inference.output_prefix=outputs/design \\
  inference.num_designs=100 \\
  potentials.guidance_scale=3.0
</code></pre>
<h3>4.2 ProteinMPNN（2022）</h3>
<p><strong>序列设计神经网络</strong>：给定骨架 → 序列，与 RFdiffusion 形成"生成-编码"闭环：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ProteinMPNN 序列设计</span>
python protein_mpnn_run.py \\
  --pdbpath scaffold.pdb \\
  --out_folder mpnn_out \\
  --num_seq_per_target 8 \\
  --batch_size 4
</code></pre>
<p><strong>标准工作流（Baker Lab 范式）</strong>：</p>
<pre><code>① RFdiffusion：生成骨架
   ↓
② ProteinMPNN：骨架 → 多组候选序列（采样温度 0.1-0.3）
   ↓
③ AF2 反向验证：序列 → 预测结构 → 与骨架 RMSD 比对
   ↓
④ 筛选高 pLDDT / 低 RMSD 候选
   ↓
⑤ 实验表达验证（酵母/噬菌体展示、ELISA）
</code></pre>
<h3>4.3 其他生成式工具</h3>
<ul>
<li><strong>Chroma</strong>：另一个扩散模型（生成 + 序列设计一体）</li>
<li><strong>Genie</strong>：更快的蛋白质扩散模型</li>
<li><strong>ESMFold</strong>（Meta，2022）：无 MSA 的极速预测（每序列 &#x3C; 1 秒），适合大规模设计筛选</li>
<li><strong>ESM3</strong>：多模态生成（序列+结构+功能）</li>
</ul>
<h2>5. binder 设计专题</h2>
<p>Binder（结合蛋白）设计是当前最活跃的应用场景：</p>
<h3>5.1 经典流程（结合现有工具）</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 靶点表面热点分析（Rosetta）</span>
<span class="hljs-comment"># 2. 热区引导的 binder 骨架生成（RFdiffusion / 或 hotspot 约束）</span>
<span class="hljs-comment"># 3. ProteinMPNN 序列设计</span>
<span class="hljs-comment"># 4. AF2 复合物预测验证（binder + target 对接）</span>
<span class="hljs-comment"># 5. 实验筛选（酵母展示、FACS 分选）</span>
</code></pre>
<h3>5.2 关键考量</h3>
<ul>
<li>结合热点残基（hotspot）的确定：实验突变扫描（alanine scan）或计算（Rosetta ddG）</li>
<li>结合亲和力预测：<code>FoldX</code>、<code>Rosetta interface_ddg</code>、AF2 PAE</li>
<li>多轮"设计-验证"迭代：实验反馈回到计算</li>
</ul>
<h2>6. 工具选择速查表</h2>
<table>
<thead>
<tr>
<th>任务</th>
<th>首选工具</th>
<th>备选</th>
</tr>
</thead>
<tbody>
<tr>
<td>结构预测</td>
<td>AlphaFold2/3（ColabFold）</td>
<td>ESMFold（速度）</td>
</tr>
<tr>
<td>骨架生成</td>
<td>RFdiffusion</td>
<td>Chroma</td>
</tr>
<tr>
<td>序列设计</td>
<td>ProteinMPNN</td>
<td>Rosetta FixBB</td>
</tr>
<tr>
<td>热稳定性设计</td>
<td>Rosetta（ddg）</td>
<td>FoldX</td>
</tr>
<tr>
<td>结合亲和力评估</td>
<td>Rosetta interface_ddg</td>
<td>FoldX / AF2 PAE</td>
</tr>
<tr>
<td>界面设计</td>
<td>Rosetta 协议</td>
<td>RFdiffusion + MPNN</td>
</tr>
<tr>
<td>大规模筛选</td>
<td>ESMFold + ProteinMPNN</td>
<td>—</td>
</tr>
</tbody>
</table>
<h2>7. 实验验证闭环</h2>
<p>计算设计必须回到实验：</p>
<pre><code>表达（E. coli / 酵母）→ 纯化 → 表征
  → SEC / DSC（稳定性）
  → 表面等离子共振 SPR / BLI（亲和力）
  → 结构验证（晶体 / 冷冻电镜 / AlphaFold 对照）
</code></pre>
<p><strong>设计成功率参考</strong>：文献报道 AF2 引导的 RFdiffusion binder 设计，实验阳性率可达 10–20%（远超传统方法），但每个成功案例背后都有多轮迭代。</p>
<h2>8. 小结</h2>
<ul>
<li><strong>经典范式</strong>：Rosetta 能量函数（可解释、可调控）</li>
<li><strong>预测范式</strong>：AlphaFold 系列（验证、复合物预测）</li>
<li><strong>生成范式</strong>：RFdiffusion + ProteinMPNN + AF2 验证闭环（当前主流）</li>
<li>设计成功的衡量标准始终是<strong>实验验证</strong></li>
</ul>
<p>下一篇将介绍虚拟筛选的主流工具：把小分子与靶点对接的完整工具链。</p>`,f_=`<h1>Mainstream Tools for Virtual Screening</h1>
<p>Virtual screening (VS) uses computation to search for potential active molecules in million-scale compound libraries and is a core approach in the early stage of drug discovery. This article introduces mainstream docking software, the complete workflow, and important considerations.</p>
<h2>1. Basic Logic of Virtual Screening</h2>
<pre><code>Receptor structure (protein/nucleic acid)
   ↓
① Receptor preparation (add hydrogens, protonation, grid)
   ↓
② Compound library (ZINC / Enamine / self-built library)
   ↓
③ Ligand preparation (3D conformation, protonation, force field parameters)
   ↓
④ Molecular docking (conformational sampling + scoring)
   ↓
⑤ Ranking and selection (Top 100-1000)
   ↓
⑥ Experimental validation (IC50 / activity assay)
</code></pre>
<h2>2. Mainstream Docking Software</h2>
<h3>2.1 AutoDock Vina (Top Open-Source Choice)</h3>
<p>The most popular open-source docking software, offering a good balance between speed and accuracy:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Receptor and ligand preparation (using MGLTools or scripts)</span>
<span class="hljs-comment">#    Receptor: pdb → pdbqt (add hydrogens, merge nonpolar hydrogens)</span>
<span class="hljs-comment">#    Ligand: sdf/mol2 → pdbqt</span>

<span class="hljs-comment"># 2. Configuration file (conf.txt)</span>
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 10.5
center_y = 20.3
center_z = -5.8
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 16
num_modes = 9

<span class="hljs-comment"># 3. Run</span>
vina --config conf.txt --out results.pdbqt --<span class="hljs-built_in">log</span> log.txt

<span class="hljs-comment"># 4. Interpretation: the lower the affinity (kcal/mol), the better</span>
<span class="hljs-comment">#  mode |   affinity | dist from best mode</span>
<span class="hljs-comment"># -----+------------+----------------------</span>
<span class="hljs-comment">#    1 |   -9.2     | 0.000</span>
<span class="hljs-comment">#    2 |   -8.7     | 2.315</span>
</code></pre>
<h3>2.2 Glide (Schrödinger, Commercial)</h3>
<p>The standard in the pharmaceutical industry:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Prepare workflow (Maestro GUI or command line)</span>
glide -overwrite -adjust \\
  -receptor receptor.maegz \\
  -ligand ligands.maegz \\
  -p glide_docking.inp \\
  -NJOBS 16

<span class="hljs-comment"># Three precision modes</span>
<span class="hljs-comment"># HTVS (high-throughput, fast) → SP (standard precision) → XP (high precision, slow)</span>
</code></pre>
<p><strong>Advantages</strong>: Accurate flexible ligand handling, well-calibrated scoring, and support from the large pharma ecosystem.</p>
<h3>2.3 DOCK 6 (Open Source)</h3>
<p>Classic geometric matching algorithm:</p>
<pre><code class="hljs language-bash">dock6 -i dock.in -o dock.out
<span class="hljs-comment"># Sphere generation → orientation search → scoring</span>
<span class="hljs-comment"># Strengths: binding pocket shape matching, initial screening of large libraries</span>
</code></pre>
<h3>2.4 Other Commonly Used Software</h3>
<table>
<thead>
<tr>
<th>Software</th>
<th>Features</th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>AutoDock4</strong></td>
<td>Classic genetic algorithm, high accuracy but slow</td>
</tr>
<tr>
<td><strong>rDock</strong></td>
<td>Open source, suitable for large libraries, easy to script</td>
</tr>
<tr>
<td><strong>GOLD</strong></td>
<td>Genetic algorithm, developed by CCDC</td>
</tr>
<tr>
<td><strong>LeDock</strong></td>
<td>Fast and easy to use (free for academia)</td>
</tr>
<tr>
<td><strong>Plants</strong></td>
<td>Based on Ant Colony Optimization</td>
</tr>
<tr>
<td><strong>smina</strong></td>
<td>Enhanced version of Vina (supports custom scoring)</td>
</tr>
</tbody>
</table>
<h2>3. Receptor Preparation (Critical First Step)</h2>
<h3>3.1 Structure Sources</h3>
<ul>
<li>Crystal/cryo-EM structures (PDB)</li>
<li>AlphaFold-predicted structures (when no experimental structure is available)</li>
<li>Note: handle flexibility (side chains) near the binding pocket</li>
</ul>
<h3>3.2 Key Preparation Points</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Common tools</span>
<span class="hljs-comment"># Schrödinger Protein Prep Wizard (commercial)</span>
<span class="hljs-comment"># ADFRsuite / MGLTools prepare_receptor4.py (free)</span>

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

<span class="hljs-comment"># Key steps</span>
<span class="hljs-comment"># 1. Remove water molecules (except conserved waters)</span>
<span class="hljs-comment"># 2. Add hydrogens, determine protonation states (around pH 7.4)</span>
<span class="hljs-comment"># 3. Assign bond orders and charges (Amber/CHARMM force fields)</span>
<span class="hljs-comment"># 4. Define the binding pocket (known active site or cavity detection: fpocket / DoGSiteScorer)</span>
</code></pre>
<h2>4. Ligand Libraries and Preparation</h2>
<h3>4.1 Compound Library Sources</h3>
<table>
<thead>
<tr>
<th>Library</th>
<th>Scale</th>
<th>Features</th>
</tr>
</thead>
<tbody>
<tr>
<td>ZINC20</td>
<td>Billions</td>
<td>Free, downloadable, with 3D conformations</td>
</tr>
<tr>
<td>Enamine REAL</td>
<td>> 40 billion</td>
<td>Highly synthesizable, commercially available</td>
</tr>
<tr>
<td>ChEMBL</td>
<td>Bioactivity data</td>
<td>Known active compounds</td>
</tr>
<tr>
<td>Self-built library</td>
<td>Custom</td>
<td>Derivative design</td>
</tr>
</tbody>
</table>
<h3>4.2 Ligand Preparation</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 3D conformation generation</span>
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

<span class="hljs-comment"># Multiple conformations (flexibility is important, especially for ring systems)</span>
obabel ligand.sdf -O conf.sdf --conformer --nconf 50

<span class="hljs-comment"># Protonation and charges</span>
<span class="hljs-comment"># Schrödinger LigPrep / Open Babel / RDKit</span>

<span class="hljs-comment"># RDKit preparation (Python)</span>
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

mol = Chem.MolFromSmiles(<span class="hljs-string">"CC(=O)Oc1ccccc1C(=O)O"</span>)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
</code></pre>
<h2>5. Scoring Functions</h2>
<h3>5.1 Three Types of Scoring Functions</h3>
<table>
<thead>
<tr>
<th>Type</th>
<th>Representatives</th>
<th>Principle</th>
</tr>
</thead>
<tbody>
<tr>
<td>Force-field-based</td>
<td>DOCK, AutoDock4</td>
<td>van der Waals + electrostatics + internal energy</td>
</tr>
<tr>
<td>Empirical</td>
<td>Glide, Vina</td>
<td>Statistical fitting to experimental data</td>
</tr>
<tr>
<td>Knowledge-based</td>
<td>GOLD/ASP</td>
<td>Statistics of atom-pair distance distributions</td>
</tr>
</tbody>
</table>
<h3>5.2 Limitations of Scoring Functions (Must-Know)</h3>
<ul>
<li><strong>Not precise</strong>: correlation coefficients are typically 0.3–0.6; they can only rank, not provide absolute affinities</li>
<li>Entropy effects, solvent effects, and induced fit are difficult to describe accurately</li>
<li><strong>Recommendation</strong>: docking scores + cross-validation with multiple tools + refined calculations using molecular dynamics (MM/PBSA, FEP)</li>
</ul>
<h2>6. Complete High-Throughput Screening Workflow (Hands-On HTS)</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Library preparation: split SDF + convert to pdbqt</span>
<span class="hljs-built_in">mkdir</span> -p ligands_pdbqt
obabel library.sdf -O ligands_pdbqt/lig_.pdbqt -m -p 7.4

<span class="hljs-comment"># 2. Batch docking (smina/Vina multi-file mode)</span>
smina -r receptor.pdbqt -l library.sdf \\
  --center_x 10.5 --center_y 20.3 --center_z -5.8 \\
  --size_x 22 --size_y 22 --size_z 22 \\
  --num_modes 1 --cpu 32 \\
  -o docked.sdf --<span class="hljs-built_in">log</span> scores.txt

<span class="hljs-comment"># 3. Ranking (by affinity)</span>
<span class="hljs-built_in">sort</span> -k2 -n scores.txt | <span class="hljs-built_in">head</span> -100
</code></pre>
<h3>6.1 Multi-Stage Funnel Strategy</h3>
<pre><code>Coarse screening (HTVS / Vina, million-scale) → top 5-10%
  ↓
Standard precision (SP / multiple conformations) → top 10-20%
  ↓
High precision (XP / induced fit) → Top 100-1000
  ↓
Consensus scoring (cross-validation with multiple software tools) → Top 50-200
  ↓
Visual inspection (reasonableness of binding modes) → Top 10-50
  ↓
Experimental testing
</code></pre>
<h3>6.2 Consensus Scoring</h3>
<p>Scoring functions from different software tools complement each other; taking intersections/ranking averages can significantly improve enrichment rates:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Example: take the intersection of Top 100 from Vina + Glide + LeDock</span>
<span class="hljs-comment"># Or average after normalizing ranks</span>
</code></pre>
<h2>7. Validation and Advanced Methods</h2>
<h3>7.1 Benchmarking</h3>
<ul>
<li><strong>DUD-E</strong>: a standard dataset (containing actives and decoys) used to evaluate a software's ability to discriminate</li>
<li>Metrics: AUC, enrichment factors (EF1%, EF5%)</li>
</ul>
<h3>7.2 Molecular Dynamics Post-Processing</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Docking results → MD simulation (GROMACS/AMBER) → MM/PBSA binding free energy</span>
gmx mdrun -deffnm complex
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ct traj.xtc ...
</code></pre>
<h3>7.3 Combining Structure-Based and Ligand-Based Approaches</h3>
<ul>
<li><strong>Structure-based</strong> (SBDD): docking + MD + FEP</li>
<li><strong>Ligand-based</strong> (LBDD): pharmacophore, QSAR, similarity search</li>
<li>Combination strategies are the mainstream of modern virtual screening</li>
</ul>
<h2>8. Common Pitfalls</h2>
<table>
<thead>
<tr>
<th>Pitfall</th>
<th>Consequence</th>
<th>Avoidance</th>
</tr>
</thead>
<tbody>
<tr>
<td>Incorrect pocket definition</td>
<td>Entire library misaligned</td>
<td>Use a known active site or experimental information</td>
</tr>
<tr>
<td>Rigid receptor</td>
<td>Misses induced fit</td>
<td>Flexible side chains / ensemble docking</td>
</tr>
<tr>
<td>Incorrect protonation states</td>
<td>Distorted electrostatic scoring</td>
<td>pH-dependent preparation (pH 7.4)</td>
</tr>
<tr>
<td>Overconfidence in scoring functions</td>
<td>High false-positive rate</td>
<td>Consensus scoring + visual inspection</td>
</tr>
<tr>
<td>Ignoring synthesizability</td>
<td>Hits cannot be synthesized</td>
<td>Filter PAINS and drug-likeness rules (Lipinski/Veber)</td>
</tr>
</tbody>
</table>
<h2>9. Summary</h2>
<ul>
<li>Three tiers: <strong>Vina</strong> (open source, fast) → <strong>Glide</strong> (commercial standard) → <strong>DOCK/rDock</strong> (large libraries)</li>
<li>Workflow: receptor preparation → ligand library → docking → multi-stage funnel → experimental validation</li>
<li>Scoring functions can only rank; consensus scoring + MD refinement improve reliability</li>
<li>The end point of virtual screening is always experimentation</li>
</ul>
<p>This completes the computational biology direction: protein design tools + virtual screening tools, with two main threads forming a closed loop of "design-screen" computational drug discovery.</p>`,g_=`<h1>虚拟筛选主流工具</h1>
<p>虚拟筛选（Virtual Screening, VS）通过计算在百万级化合物库中寻找潜在活性分子，是药物发现早期阶段的核心手段。本文介绍主流对接软件、完整流程与注意事项。</p>
<h2>1. 虚拟筛选的基本逻辑</h2>
<pre><code>受体结构（蛋白/核酸）
   ↓
① 受体准备（加氢、质子化、网格）
   ↓
② 化合物库（ZINC / Enamine / 自建库）
   ↓
③ 配体准备（3D 构象、质子化、力场参数）
   ↓
④ 分子对接（构象采样 + 打分）
   ↓
⑤ 排序筛选（Top 100-1000）
   ↓
⑥ 实验验证（IC50 / 活性测定）
</code></pre>
<h2>2. 主流对接软件</h2>
<h3>2.1 AutoDock Vina（开源首选）</h3>
<p>最流行的开源对接软件，速度与精度平衡良好：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 受体与配体准备（用 MGLTools 或脚本）</span>
<span class="hljs-comment">#    受体：pdb → pdbqt（加氢、合并非极性氢）</span>
<span class="hljs-comment">#    配体：sdf/mol2 → pdbqt</span>

<span class="hljs-comment"># 2. 配置文件（conf.txt）</span>
receptor = receptor.pdbqt
ligand = ligand.pdbqt
center_x = 10.5
center_y = 20.3
center_z = -5.8
size_x = 22
size_y = 22
size_z = 22
exhaustiveness = 16
num_modes = 9

<span class="hljs-comment"># 3. 运行</span>
vina --config conf.txt --out results.pdbqt --<span class="hljs-built_in">log</span> log.txt

<span class="hljs-comment"># 4. 解读：affinity（kcal/mol）越低越好</span>
<span class="hljs-comment">#  mode |   affinity | dist from best mode</span>
<span class="hljs-comment"># -----+------------+----------------------</span>
<span class="hljs-comment">#    1 |   -9.2     | 0.000</span>
<span class="hljs-comment">#    2 |   -8.7     | 2.315</span>
</code></pre>
<h3>2.2 Glide（Schrödinger，商业）</h3>
<p>制药行业标准：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 准备工作流（Maestro GUI 或命令行）</span>
glide -overwrite -adjust \\
  -receptor receptor.maegz \\
  -ligand ligands.maegz \\
  -p glide_docking.inp \\
  -NJOBS 16

<span class="hljs-comment"># 三种精度模式</span>
<span class="hljs-comment"># HTVS（高通量，快）→ SP（标准精度）→ XP（高精度，慢）</span>
</code></pre>
<p><strong>优势</strong>：精确的配体柔性处理、完善的打分校准、大型药企生态支持。</p>
<h3>2.3 DOCK 6（开源）</h3>
<p>经典几何匹配算法（shape matching）：</p>
<pre><code class="hljs language-bash">dock6 -i dock.in -o dock.out
<span class="hljs-comment"># 球集（sphere）生成 → 取向搜索 → 打分</span>
<span class="hljs-comment"># 擅长：结合口袋形状匹配、大库初筛</span>
</code></pre>
<h3>2.4 其他常用软件</h3>
<table>
<thead>
<tr>
<th>软件</th>
<th>特点</th>
</tr>
</thead>
<tbody>
<tr>
<td><strong>AutoDock4</strong></td>
<td>经典遗传算法，精度高速度慢</td>
</tr>
<tr>
<td><strong>rDock</strong></td>
<td>开源、适合大库、易脚本化</td>
</tr>
<tr>
<td><strong>GOLD</strong></td>
<td>遗传算法，CCDC 出品</td>
</tr>
<tr>
<td><strong>LeDock</strong></td>
<td>快速、易用（学术免费）</td>
</tr>
<tr>
<td><strong>Plants</strong></td>
<td>基于 Ant Colony 优化</td>
</tr>
<tr>
<td><strong>smina</strong></td>
<td>Vina 的增强版（支持自定义打分）</td>
</tr>
</tbody>
</table>
<h2>3. 受体准备（关键第一步）</h2>
<h3>3.1 结构来源</h3>
<ul>
<li>晶体/冷冻电镜结构（PDB）</li>
<li>AlphaFold 预测结构（无实验结构时）</li>
<li>注意：结合口袋附近的柔性（侧链）处理</li>
</ul>
<h3>3.2 准备要点</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 常用工具</span>
<span class="hljs-comment"># Schrödinger Protein Prep Wizard（商业）</span>
<span class="hljs-comment"># ADFRsuite / MGLTools prepare_receptor4.py（免费）</span>

python prepare_receptor4.py -r receptor.pdb -o receptor.pdbqt

<span class="hljs-comment"># 关键步骤</span>
<span class="hljs-comment"># 1. 去除水分子（保守水除外）</span>
<span class="hljs-comment"># 2. 加氢、确定质子化状态（pH 7.4 附近）</span>
<span class="hljs-comment"># 3. 分配键序与电荷（Amber/CHARMM 力场）</span>
<span class="hljs-comment"># 4. 定义结合口袋（已知活性位点 or 空腔检测：fpocket / DoGSiteScorer）</span>
</code></pre>
<h2>4. 配体库与准备</h2>
<h3>4.1 化合物库来源</h3>
<table>
<thead>
<tr>
<th>库</th>
<th>规模</th>
<th>特点</th>
</tr>
</thead>
<tbody>
<tr>
<td>ZINC20</td>
<td>数十亿</td>
<td>免费、可下载、带 3D 构象</td>
</tr>
<tr>
<td>Enamine REAL</td>
<td>> 400 亿</td>
<td>可合成性高、商购</td>
</tr>
<tr>
<td>ChEMBL</td>
<td>活性数据</td>
<td>已知活性化合物</td>
</tr>
<tr>
<td>自建库</td>
<td>定制</td>
<td>衍生物设计</td>
</tr>
</tbody>
</table>
<h3>4.2 配体准备</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 3D 构象生成</span>
obabel ligand.smi -O ligand.sdf --gen3d -p 7.4

<span class="hljs-comment"># 多构象（柔性重要，尤其环状体系）</span>
obabel ligand.sdf -O conf.sdf --conformer --nconf 50

<span class="hljs-comment"># 质子化与电荷</span>
<span class="hljs-comment"># Schrödinger LigPrep / Open Babel / RDKit</span>

<span class="hljs-comment"># RDKit 准备（Python）</span>
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

mol = Chem.MolFromSmiles(<span class="hljs-string">"CC(=O)Oc1ccccc1C(=O)O"</span>)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDG())
</code></pre>
<h2>5. 打分函数（Scoring Functions）</h2>
<h3>5.1 三类打分函数</h3>
<table>
<thead>
<tr>
<th>类型</th>
<th>代表</th>
<th>原理</th>
</tr>
</thead>
<tbody>
<tr>
<td>力场型</td>
<td>DOCK、AutoDock4</td>
<td>范德华+静电+内能</td>
</tr>
<tr>
<td>经验型</td>
<td>Glide、Vina</td>
<td>统计拟合实验数据</td>
</tr>
<tr>
<td>知识型</td>
<td>GOLD/ASP</td>
<td>原子对距离分布统计</td>
</tr>
</tbody>
</table>
<h3>5.2 打分函数的局限（必须知道）</h3>
<ul>
<li><strong>不精确</strong>：相关系数通常 0.3–0.6，只能排序不能给绝对亲和力</li>
<li>熵效应、溶剂效应、诱导拟合难以准确描述</li>
<li><strong>建议</strong>：对接分数 + 多工具交叉验证 + 分子动力学（MM/PBSA、FEP）精算</li>
</ul>
<h2>6. 完整高通量筛选流程（HTS 实操）</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 库准备：SDF 拆分 + 转 pdbqt</span>
<span class="hljs-built_in">mkdir</span> -p ligands_pdbqt
obabel library.sdf -O ligands_pdbqt/lig_.pdbqt -m -p 7.4

<span class="hljs-comment"># 2. 批量对接（smina/Vina 多文件模式）</span>
smina -r receptor.pdbqt -l library.sdf \\
  --center_x 10.5 --center_y 20.3 --center_z -5.8 \\
  --size_x 22 --size_y 22 --size_z 22 \\
  --num_modes 1 --cpu 32 \\
  -o docked.sdf --<span class="hljs-built_in">log</span> scores.txt

<span class="hljs-comment"># 3. 排序（按 affinity）</span>
<span class="hljs-built_in">sort</span> -k2 -n scores.txt | <span class="hljs-built_in">head</span> -100
</code></pre>
<h3>6.1 多阶段漏斗策略</h3>
<pre><code>粗筛（HTVS / Vina，百万级）→ 前 5-10%
  ↓
标准精度（SP / 多构象）→ 前 10-20%
  ↓
高精度（XP / 诱导拟合）→ Top 100-1000
  ↓
共识打分（多个软件交叉）→ Top 50-200
  ↓
视觉检查（结合模式合理性）→ Top 10-50
  ↓
实验测试
</code></pre>
<h3>6.2 共识打分（Consensus Scoring）</h3>
<p>不同软件打分函数互补，取交集/排名平均可显著提升富集率：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 示例：Vina + Glide + LeDock 三软件 Top 100 取交集</span>
<span class="hljs-comment"># 或按排名归一化后取平均</span>
</code></pre>
<h2>7. 验证与高级方法</h2>
<h3>7.1 基准测试（Benchmark）</h3>
<ul>
<li><strong>DUD-E</strong>：标准数据集（含活性物与诱饵 decoys），评估软件区分能力</li>
<li>指标：AUC、富集因子（EF1%、EF5%）</li>
</ul>
<h3>7.2 分子动力学后处理</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 对接结果 → MD 模拟（GROMACS/AMBER）→ MM/PBSA 结合自由能</span>
gmx mdrun -deffnm complex
gmx_MMPBSA -O -i mmpbsa.in -cs complex.tpr -ct traj.xtc ...
</code></pre>
<h3>7.3 结构基与配体基结合</h3>
<ul>
<li><strong>基于结构</strong>（SBDD）：对接 + MD + FEP</li>
<li><strong>基于配体</strong>（LBDD）：药效团、QSAR、相似性搜索</li>
<li>组合策略是现代虚拟筛选主流</li>
</ul>
<h2>8. 常见陷阱</h2>
<table>
<thead>
<tr>
<th>陷阱</th>
<th>后果</th>
<th>规避</th>
</tr>
</thead>
<tbody>
<tr>
<td>口袋定义错误</td>
<td>全库错位</td>
<td>用已知活性位点或实验信息</td>
</tr>
<tr>
<td>受体刚性</td>
<td>错过诱导拟合</td>
<td>柔性侧链 / 集合对接（ensemble docking）</td>
</tr>
<tr>
<td>质子化状态错误</td>
<td>静电打分失真</td>
<td>pH 相关准备（pH 7.4）</td>
</tr>
<tr>
<td>打分函数过信</td>
<td>假阳性率高</td>
<td>共识打分 + 视觉检查</td>
</tr>
<tr>
<td>忽略可合成性</td>
<td>筛出无法合成</td>
<td>过滤 PAINS、类药性规则（Lipinski/Veber）</td>
</tr>
</tbody>
</table>
<h2>9. 小结</h2>
<ul>
<li>三梯队：<strong>Vina</strong>（开源快速）→ <strong>Glide</strong>（商业标准）→ <strong>DOCK/rDock</strong>（大库）</li>
<li>流程：受体准备 → 配体库 → 对接 → 多阶段漏斗 → 实验验证</li>
<li>打分函数只能排序；共识打分 + MD 精算提升可信度</li>
<li>虚拟筛选的终点永远是实验</li>
</ul>
<p>至此计算生物学方向完成：蛋白质设计工具 + 虚拟筛选工具，两条主线构成"设计-筛选"的计算药物发现闭环。</p>`,b_=`<h1>Bash Programming: Variables, Conditionals, Loops, and Practical Scripts</h1>
<p>Bash scripts organize command-line operations into reusable programs and are a core tool for automated batch processing. This article covers variables, conditionals, loops, functions, and common tasks, from script basics to practical applications.</p>
<h2>1. Script Basics</h2>
<h3>1.1 Script Structure</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># The first line shebang: specifies the interpreter</span>

<span class="hljs-comment"># Give the script execute permission</span>
<span class="hljs-built_in">chmod</span> +x myscript.sh
<span class="hljs-comment"># Run</span>
./myscript.sh
<span class="hljs-comment"># Or run directly with bash (no execute permission needed)</span>
bash myscript.sh
</code></pre>
<h3>1.2 Debugging</h3>
<pre><code class="hljs language-bash">bash -x script.sh    <span class="hljs-comment"># Trace execution, display each command</span>
<span class="hljs-comment"># Enable inside the script: set -x (trace) / set -e (exit on error) / set -u (error on undefined variable)</span>
</code></pre>
<p>For production scripts, it's recommended to add at the beginning:</p>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">set</span> -euo pipefail   <span class="hljs-comment"># Exit on error, error on undefined variables, pipe failure detection</span>
</code></pre>
<h2>2. Variables</h2>
<h3>2.1 Definition and Reference</h3>
<pre><code class="hljs language-bash">name=<span class="hljs-string">"zorrooz"</span>        <span class="hljs-comment"># Note: no spaces around =</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>"</span>   <span class="hljs-comment"># Variable expansion inside double quotes</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">'Hello, $name'</span>   <span class="hljs-comment"># Literal output inside single quotes</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Length: <span class="hljs-variable">\${#name}</span>"</span> <span class="hljs-comment"># String length</span>
</code></pre>
<h3>2.2 Command Substitution and Arithmetic</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Command substitution: $(...)</span>
today=$(<span class="hljs-built_in">date</span> +%F)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Today is <span class="hljs-variable">$today</span>"</span>
files=$(<span class="hljs-built_in">ls</span> | <span class="hljs-built_in">wc</span> -l)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Total: <span class="hljs-variable">$files</span> files"</span>

<span class="hljs-comment"># Arithmetic: $((...))</span>
a=10
b=3
<span class="hljs-built_in">echo</span> $((a + b))    <span class="hljs-comment"># 13</span>
<span class="hljs-built_in">echo</span> $((a % b))    <span class="hljs-comment"># 1</span>
</code></pre>
<h3>2.3 Positional Parameters</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Script name: <span class="hljs-variable">$0</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"First argument: <span class="hljs-variable">$1</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Second argument: <span class="hljs-variable">$2</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"All arguments: <span class="hljs-variable">$@</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Argument count: <span class="hljs-variable">$#</span>"</span>

<span class="hljs-comment"># Run: ./script.sh file1 file2</span>
</code></pre>
<h3>2.4 Special Variables</h3>
<pre><code class="hljs language-bash">$?    <span class="hljs-comment"># Exit code of the last command (0 success, non-zero failure)</span>
$$    <span class="hljs-comment"># PID of the current process</span>
$!    <span class="hljs-comment"># PID of the background process</span>
</code></pre>
<h2>3. Conditional Statements</h2>
<h3>3.1 if / elif / else</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
score=<span class="hljs-variable">$1</span>

<span class="hljs-keyword">if</span> (( score >= <span class="hljs-number">90</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Excellent"</span>
<span class="hljs-keyword">elif</span> (( score >= <span class="hljs-number">60</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Pass"</span>
<span class="hljs-keyword">else</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Fail"</span>
<span class="hljs-keyword">fi</span>
</code></pre>
<h3>3.2 test Expressions</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># File tests</span>
[ -f file.txt ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Is a regular file"</span>
[ -d <span class="hljs-built_in">dir</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Is a directory"</span>
[ -e path ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Exists"</span>
[ -s file ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Is a non-empty file"</span>

<span class="hljs-comment"># String comparison</span>
[ <span class="hljs-string">"<span class="hljs-variable">$name</span>"</span> = <span class="hljs-string">"zorrooz"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Match"</span>
[ -z <span class="hljs-string">"<span class="hljs-variable">$var</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Variable is empty"</span>     <span class="hljs-comment"># -n non-empty</span>

<span class="hljs-comment"># Numeric comparison (inside [ ])</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"a > b"</span>    <span class="hljs-comment"># -eq -ne -gt -ge -lt -le</span>

<span class="hljs-comment"># Logical combination</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; [ -s <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"Non-empty file"</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">echo</span> <span class="hljs-string">"File does not exist"</span>
</code></pre>
<h3>3.3 case Branches</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-keyword">case</span> <span class="hljs-string">"<span class="hljs-variable">$1</span>"</span> <span class="hljs-keyword">in</span>
  start)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Start the service"</span>
    ;;
  stop|<span class="hljs-built_in">kill</span>)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Stop the service"</span>
    ;;
  *)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Usage: <span class="hljs-variable">$0</span> {start|stop}"</span>
    <span class="hljs-built_in">exit</span> 1
    ;;
<span class="hljs-keyword">esac</span>
</code></pre>
<h2>4. Loops</h2>
<h3>4.1 for Loops</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Iterate over a list</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> apple banana cherry; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Processing <span class="hljs-variable">$fruit</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># Iterate over a range</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..5}; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Iteration <span class="hljs-variable">$i</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># Iterate over files (most common)</span>
<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.fasta; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Processing <span class="hljs-variable">$file</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># C-style</span>
<span class="hljs-keyword">for</span> ((i = <span class="hljs-number">0</span>; i &#x3C; <span class="hljs-number">10</span>; i++)); <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.2 while Loops</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Read a file line by line</span>
<span class="hljs-keyword">while</span> IFS= <span class="hljs-built_in">read</span> -r line; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Read: <span class="hljs-variable">$line</span>"</span>
<span class="hljs-keyword">done</span> &#x3C; input.txt

<span class="hljs-comment"># Counting loop</span>
count=0
<span class="hljs-keyword">while</span> [ <span class="hljs-variable">$count</span> -lt 5 ]; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$count</span>
    count=$((count + <span class="hljs-number">1</span>))
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.3 break / continue</h3>
<pre><code class="hljs language-bash"><span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..10}; <span class="hljs-keyword">do</span>
    [ <span class="hljs-variable">$i</span> -eq 3 ] &#x26;&#x26; <span class="hljs-built_in">continue</span>   <span class="hljs-comment"># Skip 3</span>
    [ <span class="hljs-variable">$i</span> -eq 7 ] &#x26;&#x26; <span class="hljs-built_in">break</span>      <span class="hljs-comment"># Terminate at 7</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h2>5. Functions</h2>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>

<span class="hljs-comment"># Define a function</span>
<span class="hljs-function"><span class="hljs-title">greet</span></span>() {
    <span class="hljs-built_in">local</span> name=<span class="hljs-variable">$1</span>           <span class="hljs-comment"># local: local variable</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>!"</span>
    <span class="hljs-built_in">return</span> 0                <span class="hljs-comment"># Return value (exit code)</span>
}

<span class="hljs-comment"># Call</span>
greet <span class="hljs-string">"zorrooz"</span>

<span class="hljs-comment"># Function returning a value (note: not return, but echo capture)</span>
<span class="hljs-function"><span class="hljs-title">add</span></span>() {
    <span class="hljs-built_in">echo</span> $(( <span class="hljs-variable">$1</span> + <span class="hljs-variable">$2</span> ))
}
result=$(add 3 4)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"3 + 4 = <span class="hljs-variable">$result</span>"</span>
</code></pre>
<h2>6. Arrays and Strings</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Arrays</span>
genes=(<span class="hljs-string">"TP53"</span> <span class="hljs-string">"BRCA1"</span> <span class="hljs-string">"EGFR"</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[0]}</span>"</span>          <span class="hljs-comment"># TP53</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[@]}</span>"</span>          <span class="hljs-comment"># All elements</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${#genes[@]}</span>"</span>         <span class="hljs-comment"># Length</span>
genes+=(<span class="hljs-string">"MYC"</span>)              <span class="hljs-comment"># Append</span>

<span class="hljs-comment"># String processing</span>
<span class="hljs-built_in">seq</span>=<span class="hljs-string">"ATGCCGTAA"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq:0:3}</span>"</span>           <span class="hljs-comment"># Extract the first 3 characters: ATG</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq//T/U}</span>"</span>          <span class="hljs-comment"># Replace all T→U</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq/AT/XX}</span>"</span>         <span class="hljs-comment"># Replace the first AT</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq^^}</span>"</span>             <span class="hljs-comment"># Convert to uppercase (bash 4+)</span>
</code></pre>
<h2>7. Practical Scripts</h2>
<h3>7.1 Batch Rename</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># Change all .txt extensions to .log</span>
<span class="hljs-built_in">set</span> -euo pipefail

<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.txt; <span class="hljs-keyword">do</span>
    [ -e <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">continue</span>          <span class="hljs-comment"># Skip when there is no match</span>
    <span class="hljs-built_in">mv</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${file%.txt}</span>.log"</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Renamed: <span class="hljs-variable">$file</span> -> <span class="hljs-variable">\${file%.txt}</span>.log"</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>7.2 Bioinformatics Batch Processing: Running FASTQ Quality Control</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># Usage: ./qc_pipeline.sh data_dir output_dir</span>
<span class="hljs-built_in">set</span> -euo pipefail

DATA_DIR=<span class="hljs-variable">\${1:?"Please provide the data directory"}</span>
OUT_DIR=<span class="hljs-variable">\${2:?"Please provide the output directory"}</span>
<span class="hljs-built_in">mkdir</span> -p <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span>

<span class="hljs-keyword">for</span> fq <span class="hljs-keyword">in</span> <span class="hljs-string">"<span class="hljs-variable">$DATA_DIR</span>"</span>/*_R1.fastq.gz; <span class="hljs-keyword">do</span>
    base=$(<span class="hljs-built_in">basename</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span> _R1.fastq.gz)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"=== Processing sample: <span class="hljs-variable">$base</span> ==="</span>

    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span>
    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${DATA_DIR}</span>/<span class="hljs-variable">\${base}</span>_R2.fastq.gz"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-built_in">echo</span> <span class="hljs-string">"All done, results are in <span class="hljs-variable">$OUT_DIR</span>"</span>
</code></pre>
<h3>7.3 Log Monitoring and Alerting</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># Monitor the number of ERRORs in the log, print an alert when exceeding the threshold</span>
LOG_FILE=<span class="hljs-string">"app.log"</span>
THRESHOLD=<span class="hljs-variable">\${1:-10}</span>

count=$(grep -c <span class="hljs-string">"ERROR"</span> <span class="hljs-string">"<span class="hljs-variable">$LOG_FILE</span>"</span> || <span class="hljs-literal">true</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Current ERROR count: <span class="hljs-variable">$count</span>"</span>

<span class="hljs-keyword">if</span> [ <span class="hljs-string">"<span class="hljs-variable">$count</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$THRESHOLD</span>"</span> ]; <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"⚠️ Error count exceeds the threshold <span class="hljs-variable">$THRESHOLD</span>!"</span>
    <span class="hljs-built_in">exit</span> 1
<span class="hljs-keyword">fi</span>
</code></pre>
<h2>8. Common Pitfalls to Watch Out For</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. No spaces around the = sign in variable assignment</span>
name = <span class="hljs-string">"x"</span>    <span class="hljs-comment"># Wrong! Will be treated as a command</span>
name=<span class="hljs-string">"x"</span>      <span class="hljs-comment"># Correct</span>

<span class="hljs-comment"># 2. Quoting: paths containing spaces must be quoted</span>
file=<span class="hljs-string">"my data.txt"</span>
<span class="hljs-built_in">cat</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span>   <span class="hljs-comment"># Correct</span>
<span class="hljs-built_in">cat</span> <span class="hljs-variable">$file</span>     <span class="hljs-comment"># Wrong: splits into two arguments</span>

<span class="hljs-comment"># 3. There must be spaces on both sides inside [ ]</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> = <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ]   <span class="hljs-comment"># Correct</span>
[<span class="hljs-string">"<span class="hljs-variable">$a</span>"</span>=<span class="hljs-string">"<span class="hljs-variable">$b</span>"</span>]       <span class="hljs-comment"># Wrong</span>

<span class="hljs-comment"># 4. The pitfall of set -e in pipelines: grep returns non-zero when no match is found, which will interrupt the script</span>
grep <span class="hljs-string">"x"</span> file || <span class="hljs-literal">true</span>    <span class="hljs-comment"># Explicitly tolerate</span>
</code></pre>
<h2>9. Summary</h2>
<ul>
<li>Script header: <code>#!/bin/bash</code> + <code>set -euo pipefail</code></li>
<li>Variables: <code>$name</code>, <code>\${#name}</code>, <code>$(cmd)</code>, <code>$((...))</code></li>
<li>Conditionals: <code>[ condition ]</code>, <code>(( arithmetic ))</code>, <code>case</code></li>
<li>Loops: <code>for file in *.txt</code>, <code>while read line</code></li>
<li>Functions: <code>local</code> variables, capturing return values with echo</li>
</ul>
<p>At this point, the programming track tutorials are complete: Python ×3, R ×3, Linux ×1, Bash ×1, from zero to hands-on practice.</p>`,__=`<h1>Bash 编程：变量、条件、循环与实用脚本</h1>
<p>Bash 脚本把命令行操作组织成可复用的程序，是自动化批处理的核心工具。本文从脚本基础到实战，覆盖变量、条件、循环、函数与常见任务。</p>
<h2>1. 脚本基础</h2>
<h3>1.1 脚本结构</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 第一行 shebang：指定解释器</span>

<span class="hljs-comment"># 给脚本执行权限</span>
<span class="hljs-built_in">chmod</span> +x myscript.sh
<span class="hljs-comment"># 运行</span>
./myscript.sh
<span class="hljs-comment"># 或用 bash 直接运行（无需执行权限）</span>
bash myscript.sh
</code></pre>
<h3>1.2 调试</h3>
<pre><code class="hljs language-bash">bash -x script.sh    <span class="hljs-comment"># 跟踪执行，显示每条命令</span>
<span class="hljs-comment"># 脚本内开启：set -x（跟踪）/ set -e（出错即退出）/ set -u（未定义变量报错）</span>
</code></pre>
<p>生产脚本推荐在开头加：</p>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">set</span> -euo pipefail   <span class="hljs-comment"># 出错即停、未定义变量报错、管道失败检测</span>
</code></pre>
<h2>2. 变量</h2>
<h3>2.1 定义与引用</h3>
<pre><code class="hljs language-bash">name=<span class="hljs-string">"zorrooz"</span>        <span class="hljs-comment"># 注意：= 两边不能有空格</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>"</span>   <span class="hljs-comment"># 双引号内变量展开</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">'Hello, $name'</span>   <span class="hljs-comment"># 单引号内原样输出</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"长度: <span class="hljs-variable">\${#name}</span>"</span> <span class="hljs-comment"># 字符串长度</span>
</code></pre>
<h3>2.2 命令替换与算术</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 命令替换：$(...)</span>
today=$(<span class="hljs-built_in">date</span> +%F)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"今天是 <span class="hljs-variable">$today</span>"</span>
files=$(<span class="hljs-built_in">ls</span> | <span class="hljs-built_in">wc</span> -l)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"共 <span class="hljs-variable">$files</span> 个文件"</span>

<span class="hljs-comment"># 算术：$((...))</span>
a=10
b=3
<span class="hljs-built_in">echo</span> $((a + b))    <span class="hljs-comment"># 13</span>
<span class="hljs-built_in">echo</span> $((a % b))    <span class="hljs-comment"># 1</span>
</code></pre>
<h3>2.3 位置参数</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"脚本名: <span class="hljs-variable">$0</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"第一个参数: <span class="hljs-variable">$1</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"第二个参数: <span class="hljs-variable">$2</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"所有参数: <span class="hljs-variable">$@</span>"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"参数个数: <span class="hljs-variable">$#</span>"</span>

<span class="hljs-comment"># 运行：./script.sh file1 file2</span>
</code></pre>
<h3>2.4 特殊变量</h3>
<pre><code class="hljs language-bash">$?    <span class="hljs-comment"># 上一条命令的退出码（0 成功，非 0 失败）</span>
$$    <span class="hljs-comment"># 当前进程 PID</span>
$!    <span class="hljs-comment"># 后台进程 PID</span>
</code></pre>
<h2>3. 条件判断</h2>
<h3>3.1 if / elif / else</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
score=<span class="hljs-variable">$1</span>

<span class="hljs-keyword">if</span> (( score >= <span class="hljs-number">90</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"优秀"</span>
<span class="hljs-keyword">elif</span> (( score >= <span class="hljs-number">60</span> )); <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"及格"</span>
<span class="hljs-keyword">else</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"不及格"</span>
<span class="hljs-keyword">fi</span>
</code></pre>
<h3>3.2 test 表达式</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 文件判断</span>
[ -f file.txt ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"是普通文件"</span>
[ -d <span class="hljs-built_in">dir</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"是目录"</span>
[ -e path ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"存在"</span>
[ -s file ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"非空文件"</span>

<span class="hljs-comment"># 字符串比较</span>
[ <span class="hljs-string">"<span class="hljs-variable">$name</span>"</span> = <span class="hljs-string">"zorrooz"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"匹配"</span>
[ -z <span class="hljs-string">"<span class="hljs-variable">$var</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"变量为空"</span>     <span class="hljs-comment"># -n 非空</span>

<span class="hljs-comment"># 数值比较（[ ] 内）</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"a > b"</span>    <span class="hljs-comment"># -eq -ne -gt -ge -lt -le</span>

<span class="hljs-comment"># 逻辑组合</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; [ -s <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] &#x26;&#x26; <span class="hljs-built_in">echo</span> <span class="hljs-string">"非空文件"</span>
[ -f <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">echo</span> <span class="hljs-string">"文件不存在"</span>
</code></pre>
<h3>3.3 case 分支</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-keyword">case</span> <span class="hljs-string">"<span class="hljs-variable">$1</span>"</span> <span class="hljs-keyword">in</span>
  start)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"启动服务"</span>
    ;;
  stop|<span class="hljs-built_in">kill</span>)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"停止服务"</span>
    ;;
  *)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"用法: <span class="hljs-variable">$0</span> {start|stop}"</span>
    <span class="hljs-built_in">exit</span> 1
    ;;
<span class="hljs-keyword">esac</span>
</code></pre>
<h2>4. 循环</h2>
<h3>4.1 for 循环</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 遍历列表</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> apple banana cherry; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"处理 <span class="hljs-variable">$fruit</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># 遍历范围</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..5}; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"第 <span class="hljs-variable">$i</span> 次"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># 遍历文件（最常用）</span>
<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.fasta; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"处理 <span class="hljs-variable">$file</span>"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-comment"># C 风格</span>
<span class="hljs-keyword">for</span> ((i = <span class="hljs-number">0</span>; i &#x3C; <span class="hljs-number">10</span>; i++)); <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.2 while 循环</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 逐行读取文件</span>
<span class="hljs-keyword">while</span> IFS= <span class="hljs-built_in">read</span> -r line; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"读取: <span class="hljs-variable">$line</span>"</span>
<span class="hljs-keyword">done</span> &#x3C; input.txt

<span class="hljs-comment"># 计数循环</span>
count=0
<span class="hljs-keyword">while</span> [ <span class="hljs-variable">$count</span> -lt 5 ]; <span class="hljs-keyword">do</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$count</span>
    count=$((count + <span class="hljs-number">1</span>))
<span class="hljs-keyword">done</span>
</code></pre>
<h3>4.3 break / continue</h3>
<pre><code class="hljs language-bash"><span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> {1..10}; <span class="hljs-keyword">do</span>
    [ <span class="hljs-variable">$i</span> -eq 3 ] &#x26;&#x26; <span class="hljs-built_in">continue</span>   <span class="hljs-comment"># 跳过 3</span>
    [ <span class="hljs-variable">$i</span> -eq 7 ] &#x26;&#x26; <span class="hljs-built_in">break</span>      <span class="hljs-comment"># 到 7 终止</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-variable">$i</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h2>5. 函数</h2>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>

<span class="hljs-comment"># 定义函数</span>
<span class="hljs-function"><span class="hljs-title">greet</span></span>() {
    <span class="hljs-built_in">local</span> name=<span class="hljs-variable">$1</span>           <span class="hljs-comment"># local：局部变量</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"Hello, <span class="hljs-variable">$name</span>!"</span>
    <span class="hljs-built_in">return</span> 0                <span class="hljs-comment"># 返回值（退出码）</span>
}

<span class="hljs-comment"># 调用</span>
greet <span class="hljs-string">"zorrooz"</span>

<span class="hljs-comment"># 函数返回数值（注意：不是 return，是 echo 捕获）</span>
<span class="hljs-function"><span class="hljs-title">add</span></span>() {
    <span class="hljs-built_in">echo</span> $(( <span class="hljs-variable">$1</span> + <span class="hljs-variable">$2</span> ))
}
result=$(add 3 4)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"3 + 4 = <span class="hljs-variable">$result</span>"</span>
</code></pre>
<h2>6. 数组与字符串</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 数组</span>
genes=(<span class="hljs-string">"TP53"</span> <span class="hljs-string">"BRCA1"</span> <span class="hljs-string">"EGFR"</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[0]}</span>"</span>          <span class="hljs-comment"># TP53</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${genes[@]}</span>"</span>          <span class="hljs-comment"># 所有元素</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${#genes[@]}</span>"</span>         <span class="hljs-comment"># 长度</span>
genes+=(<span class="hljs-string">"MYC"</span>)              <span class="hljs-comment"># 追加</span>

<span class="hljs-comment"># 字符串处理</span>
<span class="hljs-built_in">seq</span>=<span class="hljs-string">"ATGCCGTAA"</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq:0:3}</span>"</span>           <span class="hljs-comment"># 截取前 3 个字符：ATG</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq//T/U}</span>"</span>          <span class="hljs-comment"># 全部替换 T→U</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq/AT/XX}</span>"</span>         <span class="hljs-comment"># 替换第一个 AT</span>
<span class="hljs-built_in">echo</span> <span class="hljs-string">"<span class="hljs-variable">\${seq^^}</span>"</span>             <span class="hljs-comment"># 转大写（bash 4+）</span>
</code></pre>
<h2>7. 实战脚本</h2>
<h3>7.1 批量重命名</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 把所有 .txt 后缀改为 .log</span>
<span class="hljs-built_in">set</span> -euo pipefail

<span class="hljs-keyword">for</span> file <span class="hljs-keyword">in</span> *.txt; <span class="hljs-keyword">do</span>
    [ -e <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> ] || <span class="hljs-built_in">continue</span>          <span class="hljs-comment"># 无匹配时跳过</span>
    <span class="hljs-built_in">mv</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${file%.txt}</span>.log"</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"重命名: <span class="hljs-variable">$file</span> -> <span class="hljs-variable">\${file%.txt}</span>.log"</span>
<span class="hljs-keyword">done</span>
</code></pre>
<h3>7.2 生信批处理：批量跑 FASTQ 质控</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 用法: ./qc_pipeline.sh data_dir output_dir</span>
<span class="hljs-built_in">set</span> -euo pipefail

DATA_DIR=<span class="hljs-variable">\${1:?"请提供数据目录"}</span>
OUT_DIR=<span class="hljs-variable">\${2:?"请提供输出目录"}</span>
<span class="hljs-built_in">mkdir</span> -p <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span>

<span class="hljs-keyword">for</span> fq <span class="hljs-keyword">in</span> <span class="hljs-string">"<span class="hljs-variable">$DATA_DIR</span>"</span>/*_R1.fastq.gz; <span class="hljs-keyword">do</span>
    base=$(<span class="hljs-built_in">basename</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span> _R1.fastq.gz)
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"=== 处理样本: <span class="hljs-variable">$base</span> ==="</span>

    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">$fq</span>"</span>
    fastqc -o <span class="hljs-string">"<span class="hljs-variable">$OUT_DIR</span>"</span> <span class="hljs-string">"<span class="hljs-variable">\${DATA_DIR}</span>/<span class="hljs-variable">\${base}</span>_R2.fastq.gz"</span>
<span class="hljs-keyword">done</span>

<span class="hljs-built_in">echo</span> <span class="hljs-string">"全部完成，结果位于 <span class="hljs-variable">$OUT_DIR</span>"</span>
</code></pre>
<h3>7.3 日志监控与告警</h3>
<pre><code class="hljs language-bash"><span class="hljs-meta">#!/bin/bash</span>
<span class="hljs-comment"># 监控日志中 ERROR 数量，超过阈值打印告警</span>
LOG_FILE=<span class="hljs-string">"app.log"</span>
THRESHOLD=<span class="hljs-variable">\${1:-10}</span>

count=$(grep -c <span class="hljs-string">"ERROR"</span> <span class="hljs-string">"<span class="hljs-variable">$LOG_FILE</span>"</span> || <span class="hljs-literal">true</span>)
<span class="hljs-built_in">echo</span> <span class="hljs-string">"当前 ERROR 数量: <span class="hljs-variable">$count</span>"</span>

<span class="hljs-keyword">if</span> [ <span class="hljs-string">"<span class="hljs-variable">$count</span>"</span> -gt <span class="hljs-string">"<span class="hljs-variable">$THRESHOLD</span>"</span> ]; <span class="hljs-keyword">then</span>
    <span class="hljs-built_in">echo</span> <span class="hljs-string">"⚠️ 错误数量超过阈值 <span class="hljs-variable">$THRESHOLD</span>！"</span>
    <span class="hljs-built_in">exit</span> 1
<span class="hljs-keyword">fi</span>
</code></pre>
<h2>8. 常见坑位提醒</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 变量赋值等号两边不能有空格</span>
name = <span class="hljs-string">"x"</span>    <span class="hljs-comment"># 错误！会被当成命令</span>
name=<span class="hljs-string">"x"</span>      <span class="hljs-comment"># 正确</span>

<span class="hljs-comment"># 2. 引号：路径含空格必须加引号</span>
file=<span class="hljs-string">"my data.txt"</span>
<span class="hljs-built_in">cat</span> <span class="hljs-string">"<span class="hljs-variable">$file</span>"</span>   <span class="hljs-comment"># 正确</span>
<span class="hljs-built_in">cat</span> <span class="hljs-variable">$file</span>     <span class="hljs-comment"># 错误：拆成两个参数</span>

<span class="hljs-comment"># 3. [ ] 内部两端必须有空格</span>
[ <span class="hljs-string">"<span class="hljs-variable">$a</span>"</span> = <span class="hljs-string">"<span class="hljs-variable">$b</span>"</span> ]   <span class="hljs-comment"># 正确</span>
[<span class="hljs-string">"<span class="hljs-variable">$a</span>"</span>=<span class="hljs-string">"<span class="hljs-variable">$b</span>"</span>]       <span class="hljs-comment"># 错误</span>

<span class="hljs-comment"># 4. 管道中 set -e 的坑：grep 无匹配返回非 0 会中断</span>
grep <span class="hljs-string">"x"</span> file || <span class="hljs-literal">true</span>    <span class="hljs-comment"># 显式容忍</span>
</code></pre>
<h2>9. 小结</h2>
<ul>
<li>脚本头：<code>#!/bin/bash</code> + <code>set -euo pipefail</code></li>
<li>变量：<code>$name</code>、<code>\${#name}</code>、<code>$(cmd)</code>、<code>$((...))</code></li>
<li>判断：<code>[ 条件 ]</code>、<code>(( 算术 ))</code>、<code>case</code></li>
<li>循环：<code>for file in *.txt</code>、<code>while read line</code></li>
<li>函数：<code>local</code> 变量、echo 捕获返回值</li>
</ul>
<p>至此编程方向教程完成：Python ×3、R ×3、Linux ×1、Bash ×1，从零到实战。</p>`,y_=`<h1>Linux Command Line Basics: File System, Permissions, and Text Processing</h1>
<p>The Linux command line (Shell) is a fundamental skill for bioinformatics and scientific computing: server operations, workflow construction, and batch processing all depend on it. This article covers the most commonly used commands, and it is recommended to practice while learning using WSL or a remote server.</p>
<h2>1. Terminal and Shell</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">echo</span> <span class="hljs-variable">$SHELL</span>        <span class="hljs-comment"># View the current shell, usually /bin/bash</span>
<span class="hljs-built_in">whoami</span>             <span class="hljs-comment"># Current user</span>
<span class="hljs-built_in">pwd</span>                <span class="hljs-comment"># Current directory (print working directory)</span>
</code></pre>
<h2>2. File System Navigation</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span>                  <span class="hljs-comment"># List files</span>
<span class="hljs-built_in">ls</span> -l               <span class="hljs-comment"># Detailed information (permissions, size, time)</span>
<span class="hljs-built_in">ls</span> -a               <span class="hljs-comment"># Include hidden files (starting with .)</span>
<span class="hljs-built_in">ls</span> -lh              <span class="hljs-comment"># Human-readable sizes</span>

<span class="hljs-built_in">cd</span> /home/user       <span class="hljs-comment"># Enter directory</span>
<span class="hljs-built_in">cd</span> ..               <span class="hljs-comment"># Parent directory</span>
<span class="hljs-built_in">cd</span> ~                <span class="hljs-comment"># User home directory</span>
<span class="hljs-built_in">cd</span> -                <span class="hljs-comment"># Previous directory</span>

<span class="hljs-built_in">pwd</span>                 <span class="hljs-comment"># Display current location</span>
</code></pre>
<p><strong>Path rules</strong>: <code>/</code> is the root directory; <code>.</code> is the current directory; <code>..</code> is the parent directory; <code>~</code> is the home directory.</p>
<h2>3. File and Directory Operations</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Create</span>
<span class="hljs-built_in">touch</span> file.txt      <span class="hljs-comment"># Create an empty file</span>
<span class="hljs-built_in">mkdir</span> data          <span class="hljs-comment"># Create a directory</span>
<span class="hljs-built_in">mkdir</span> -p a/b/c      <span class="hljs-comment"># Recursively create multi-level directories</span>

<span class="hljs-comment"># Copy / Move / Delete</span>
<span class="hljs-built_in">cp</span> file.txt copy.txt
<span class="hljs-built_in">cp</span> -r data/ data_backup/    <span class="hljs-comment"># Copy directory</span>
<span class="hljs-built_in">mv</span> file.txt newname.txt     <span class="hljs-comment"># Rename/move</span>
<span class="hljs-built_in">rm</span> file.txt                 <span class="hljs-comment"># Delete file</span>
<span class="hljs-built_in">rm</span> -r data/                 <span class="hljs-comment"># Delete directory (dangerous!)</span>
<span class="hljs-built_in">rm</span> -rf data/                <span class="hljs-comment"># Force recursive delete (use with caution)</span>

<span class="hljs-comment"># View contents</span>
<span class="hljs-built_in">cat</span> file.txt        <span class="hljs-comment"># Output all contents</span>
less file.txt       <span class="hljs-comment"># Page through (q to quit, / to search)</span>
<span class="hljs-built_in">head</span> -n 5 file.txt  <span class="hljs-comment"># First 5 lines</span>
<span class="hljs-built_in">tail</span> -n 5 file.txt  <span class="hljs-comment"># Last 5 lines</span>
<span class="hljs-built_in">tail</span> -f log.txt     <span class="hljs-comment"># Follow in real time (commonly used for viewing logs)</span>
</code></pre>
<blockquote>
<p><strong><code>rm -rf</code> is a dangerous command</strong>: confirm the path before executing, or first <code>ls</code> to verify.</p>
</blockquote>
<h2>4. Wildcards and Redirection</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Wildcards</span>
<span class="hljs-built_in">ls</span> *.fasta          <span class="hljs-comment"># All .fasta files</span>
<span class="hljs-built_in">ls</span> data_??.txt      <span class="hljs-comment"># ? matches a single character</span>
<span class="hljs-built_in">ls</span> [abc]*           <span class="hljs-comment"># Starts with a/b/c</span>

<span class="hljs-comment"># Redirection</span>
<span class="hljs-built_in">ls</span> > list.txt       <span class="hljs-comment"># Overwrite</span>
<span class="hljs-built_in">ls</span> >> list.txt      <span class="hljs-comment"># Append</span>
<span class="hljs-built_in">ls</span> 2> error.log     <span class="hljs-comment"># Redirect error output</span>
<span class="hljs-built_in">ls</span> 2>/dev/null      <span class="hljs-comment"># Discard error output</span>

<span class="hljs-comment"># Pipe: output of the previous command becomes input of the next</span>
<span class="hljs-built_in">ls</span> -l | <span class="hljs-built_in">wc</span> -l                   <span class="hljs-comment"># Count number of files</span>
<span class="hljs-built_in">history</span> | grep python           <span class="hljs-comment"># Search history commands</span>
</code></pre>
<h2>5. Permission Management</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span> -l
<span class="hljs-comment"># -rw-r--r-- 1 user group 1234 Aug 4 10:00 file.txt</span>
<span class="hljs-comment"># ^permissions       ^owner ^group</span>

<span class="hljs-built_in">chmod</span> 755 script.sh    <span class="hljs-comment"># rwxr-xr-x: owner read/write/execute, others read/execute only</span>
<span class="hljs-built_in">chmod</span> +x script.sh     <span class="hljs-comment"># Add execute permission (required to run scripts)</span>
<span class="hljs-built_in">chmod</span> -w file.txt      <span class="hljs-comment"># Remove write permission</span>

<span class="hljs-built_in">chown</span> user:group file  <span class="hljs-comment"># Change owner and group (requires root)</span>
</code></pre>
<p>Permission numbers: <code>r=4</code>, <code>w=2</code>, <code>x=1</code>, and the three digits represent owner/group/others respectively.</p>
<h2>6. The Three Musketeers of Text Processing</h2>
<h3>6.1 grep: Search</h3>
<pre><code class="hljs language-bash">grep <span class="hljs-string">"TP53"</span> genes.txt          <span class="hljs-comment"># Find lines containing TP53</span>
grep -i <span class="hljs-string">"tp53"</span> genes.txt       <span class="hljs-comment"># Ignore case</span>
grep -v <span class="hljs-string">"comment"</span> file.txt     <span class="hljs-comment"># Reverse match (exclude)</span>
grep -c <span class="hljs-string">"gene"</span> file.txt        <span class="hljs-comment"># Count</span>
grep -n <span class="hljs-string">"pattern"</span> file.txt     <span class="hljs-comment"># Show line numbers</span>
grep -r <span class="hljs-string">"TODO"</span> src/            <span class="hljs-comment"># Recursively search directory</span>

<span class="hljs-comment"># Pipe combinations</span>
ps aux | grep python           <span class="hljs-comment"># Find python processes</span>
<span class="hljs-built_in">cat</span> reads.fastq | grep -c <span class="hljs-string">"^@"</span> <span class="hljs-comment"># Count number of sequences in FASTQ</span>
</code></pre>
<h3>6.2 sed: Stream Editing</h3>
<pre><code class="hljs language-bash">sed -n <span class="hljs-string">'5,10p'</span> file.txt        <span class="hljs-comment"># Print lines 5-10</span>
sed <span class="hljs-string">'s/old/new/'</span> file.txt      <span class="hljs-comment"># Replace first match in each line</span>
sed <span class="hljs-string">'s/old/new/g'</span> file.txt     <span class="hljs-comment"># Global replace</span>
sed -i <span class="hljs-string">'s/old/new/g'</span> file.txt  <span class="hljs-comment"># Modify file in place (-i in-place)</span>
sed <span class="hljs-string">'/^#/d'</span> config.conf        <span class="hljs-comment"># Delete comment lines</span>
</code></pre>
<h3>6.3 awk: Column-Based Processing</h3>
<pre><code class="hljs language-bash">awk <span class="hljs-string">'{print $1}'</span> file.txt      <span class="hljs-comment"># Print first column</span>
awk -F<span class="hljs-string">','</span> <span class="hljs-string">'{print $1, $3}'</span> data.csv   <span class="hljs-comment"># Separate by comma</span>
awk <span class="hljs-string">'$3 > 50 {print $1}'</span> data.txt     <span class="hljs-comment"># Conditional filter</span>
awk <span class="hljs-string">'{sum += $2} END {print sum}'</span> data.txt   <span class="hljs-comment"># Sum second column</span>
awk <span class="hljs-string">'NR > 1 {print}'</span> file.txt  <span class="hljs-comment"># Skip header line</span>
</code></pre>
<blockquote>
<p>awk field separator defaults to whitespace; <code>-F','</code> specifies comma.</p>
</blockquote>
<h2>7. Process Management</h2>
<pre><code class="hljs language-bash">ps aux                <span class="hljs-comment"># View all processes</span>
ps aux | grep python  <span class="hljs-comment"># Find specific processes</span>
top                   <span class="hljs-comment"># Real-time monitoring (q to quit)</span>
htop                  <span class="hljs-comment"># Enhanced monitoring (more intuitive)</span>

<span class="hljs-built_in">kill</span> PID              <span class="hljs-comment"># Terminate process</span>
<span class="hljs-built_in">kill</span> -9 PID           <span class="hljs-comment"># Force kill</span>
pkill python          <span class="hljs-comment"># Terminate by name</span>

<span class="hljs-comment"># Background execution</span>
python script.py &#x26;    <span class="hljs-comment"># Run in background</span>
<span class="hljs-built_in">nohup</span> python script.py > log.txt 2>&#x26;1 &#x26;   <span class="hljs-comment"># Run detached from terminal (common on servers)</span>
<span class="hljs-built_in">jobs</span>                  <span class="hljs-comment"># View background jobs</span>
<span class="hljs-built_in">fg</span>                    <span class="hljs-comment"># Bring to foreground</span>
</code></pre>
<blockquote>
<p>For long-running tasks on remote servers, use <code>nohup ... &#x26;</code> or <code>tmux</code> / <code>screen</code> to avoid task interruption due to terminal disconnection.</p>
</blockquote>
<h2>8. Practical Combinations</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Count files and sizes</span>
<span class="hljs-built_in">ls</span> -lh | awk <span class="hljs-string">'{print $5, $9}'</span>
<span class="hljs-built_in">du</span> -sh data/          <span class="hljs-comment"># Total size of directory</span>
<span class="hljs-built_in">df</span> -h                 <span class="hljs-comment"># Disk space</span>

<span class="hljs-comment"># Compression and decompression</span>
tar -czf archive.tar.gz <span class="hljs-built_in">dir</span>/    <span class="hljs-comment"># Pack and compress</span>
tar -xzf archive.tar.gz         <span class="hljs-comment"># Decompress</span>
zip -r archive.zip <span class="hljs-built_in">dir</span>/
unzip archive.zip

<span class="hljs-comment"># Find files</span>
find . -name <span class="hljs-string">"*.log"</span>            <span class="hljs-comment"># Find by name</span>
find . -size +100M              <span class="hljs-comment"># Find by size</span>
<span class="hljs-built_in">which</span> python                    <span class="hljs-comment"># Path of the command</span>
</code></pre>
<h2>9. WSL and Development Environment</h2>
<p>Under Windows, WSL (Windows Subsystem for Linux) is recommended:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Install (in administrator PowerShell)</span>
wsl --install

<span class="hljs-comment"># Enter the Ubuntu subsystem</span>
wsl

<span class="hljs-comment"># Configure the development environment in WSL</span>
<span class="hljs-built_in">sudo</span> apt update
<span class="hljs-built_in">sudo</span> apt install build-essential git python3 python3-pip
</code></pre>
<p>VS Code connects to WSL: after installing the Remote - WSL extension, <code>code .</code> opens directly.</p>
<h2>10. Summary</h2>
<ul>
<li>Navigation: <code>cd</code> / <code>ls</code> / <code>pwd</code>; operations: <code>cp</code> / <code>mv</code> / <code>rm</code> / <code>mkdir</code></li>
<li>The three musketeers: <code>grep</code> for searching, <code>sed</code> for replacing, <code>awk</code> for column-based processing</li>
<li>The pipe <code>|</code> and redirection <code>></code> are the core of combining commands</li>
<li>Permissions: <code>chmod 755</code> / <code>chmod +x</code>; processes: <code>ps</code> / <code>kill</code> / <code>nohup</code></li>
</ul>
<p>The next article will introduce Bash programming: organizing commands into reusable scripts.</p>`,v_=`<h1>Linux 命令行基础：文件系统、权限与文本处理</h1>
<p>Linux 命令行（Shell）是生物信息学与科学计算的基本功：服务器操作、流程搭建、批量处理都离不开它。本文覆盖最常用的命令，建议配合 WSL 或远程服务器边学边练。</p>
<h2>1. 终端与 Shell</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">echo</span> <span class="hljs-variable">$SHELL</span>        <span class="hljs-comment"># 查看当前 shell，通常为 /bin/bash</span>
<span class="hljs-built_in">whoami</span>             <span class="hljs-comment"># 当前用户</span>
<span class="hljs-built_in">pwd</span>                <span class="hljs-comment"># 当前目录（print working directory）</span>
</code></pre>
<h2>2. 文件系统导航</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span>                  <span class="hljs-comment"># 列出文件</span>
<span class="hljs-built_in">ls</span> -l               <span class="hljs-comment"># 详细信息（权限、大小、时间）</span>
<span class="hljs-built_in">ls</span> -a               <span class="hljs-comment"># 包含隐藏文件（以 . 开头）</span>
<span class="hljs-built_in">ls</span> -lh              <span class="hljs-comment"># 人类可读大小</span>

<span class="hljs-built_in">cd</span> /home/user       <span class="hljs-comment"># 进入目录</span>
<span class="hljs-built_in">cd</span> ..               <span class="hljs-comment"># 上级目录</span>
<span class="hljs-built_in">cd</span> ~                <span class="hljs-comment"># 用户主目录</span>
<span class="hljs-built_in">cd</span> -                <span class="hljs-comment"># 上一个目录</span>

<span class="hljs-built_in">pwd</span>                 <span class="hljs-comment"># 显示当前位置</span>
</code></pre>
<p><strong>路径规则</strong>：<code>/</code> 是根目录；<code>.</code> 当前目录；<code>..</code> 上级目录；<code>~</code> 主目录。</p>
<h2>3. 文件与目录操作</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 创建</span>
<span class="hljs-built_in">touch</span> file.txt      <span class="hljs-comment"># 创建空文件</span>
<span class="hljs-built_in">mkdir</span> data          <span class="hljs-comment"># 创建目录</span>
<span class="hljs-built_in">mkdir</span> -p a/b/c      <span class="hljs-comment"># 递归创建多级目录</span>

<span class="hljs-comment"># 复制 / 移动 / 删除</span>
<span class="hljs-built_in">cp</span> file.txt copy.txt
<span class="hljs-built_in">cp</span> -r data/ data_backup/    <span class="hljs-comment"># 复制目录</span>
<span class="hljs-built_in">mv</span> file.txt newname.txt     <span class="hljs-comment"># 重命名/移动</span>
<span class="hljs-built_in">rm</span> file.txt                 <span class="hljs-comment"># 删除文件</span>
<span class="hljs-built_in">rm</span> -r data/                 <span class="hljs-comment"># 删除目录（危险！）</span>
<span class="hljs-built_in">rm</span> -rf data/                <span class="hljs-comment"># 强制递归删除（慎用）</span>

<span class="hljs-comment"># 查看内容</span>
<span class="hljs-built_in">cat</span> file.txt        <span class="hljs-comment"># 输出全部</span>
less file.txt       <span class="hljs-comment"># 分页浏览（q 退出，/搜索）</span>
<span class="hljs-built_in">head</span> -n 5 file.txt  <span class="hljs-comment"># 前 5 行</span>
<span class="hljs-built_in">tail</span> -n 5 file.txt  <span class="hljs-comment"># 后 5 行</span>
<span class="hljs-built_in">tail</span> -f log.txt     <span class="hljs-comment"># 实时跟踪（查看日志常用）</span>
</code></pre>
<blockquote>
<p><strong>rm -rf 是危险命令</strong>：执行前确认路径，或先 <code>ls</code> 验证。</p>
</blockquote>
<h2>4. 通配符与重定向</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 通配符</span>
<span class="hljs-built_in">ls</span> *.fasta          <span class="hljs-comment"># 所有 .fasta 文件</span>
<span class="hljs-built_in">ls</span> data_??.txt      <span class="hljs-comment"># ? 匹配单个字符</span>
<span class="hljs-built_in">ls</span> [abc]*           <span class="hljs-comment"># 以 a/b/c 开头</span>

<span class="hljs-comment"># 重定向</span>
<span class="hljs-built_in">ls</span> > list.txt       <span class="hljs-comment"># 覆盖写入</span>
<span class="hljs-built_in">ls</span> >> list.txt      <span class="hljs-comment"># 追加写入</span>
<span class="hljs-built_in">ls</span> 2> error.log     <span class="hljs-comment"># 错误输出重定向</span>
<span class="hljs-built_in">ls</span> 2>/dev/null      <span class="hljs-comment"># 丢弃错误输出</span>

<span class="hljs-comment"># 管道：前一个命令的输出作为后一个的输入</span>
<span class="hljs-built_in">ls</span> -l | <span class="hljs-built_in">wc</span> -l                   <span class="hljs-comment"># 统计文件数</span>
<span class="hljs-built_in">history</span> | grep python           <span class="hljs-comment"># 在历史命令中搜索</span>
</code></pre>
<h2>5. 权限管理</h2>
<pre><code class="hljs language-bash"><span class="hljs-built_in">ls</span> -l
<span class="hljs-comment"># -rw-r--r-- 1 user group 1234 Aug 4 10:00 file.txt</span>
<span class="hljs-comment"># ^权限       ^属主 ^属组</span>

<span class="hljs-built_in">chmod</span> 755 script.sh    <span class="hljs-comment"># rwxr-xr-x：属主读写执行，其他只读执行</span>
<span class="hljs-built_in">chmod</span> +x script.sh     <span class="hljs-comment"># 添加执行权限（运行脚本必需）</span>
<span class="hljs-built_in">chmod</span> -w file.txt      <span class="hljs-comment"># 去除写权限</span>

<span class="hljs-built_in">chown</span> user:group file  <span class="hljs-comment"># 修改属主属组（需要 root）</span>
</code></pre>
<p>权限数字：<code>r=4</code>、<code>w=2</code>、<code>x=1</code>，三位分别表示属主/属组/其他。</p>
<h2>6. 文本处理三剑客</h2>
<h3>6.1 grep：搜索</h3>
<pre><code class="hljs language-bash">grep <span class="hljs-string">"TP53"</span> genes.txt          <span class="hljs-comment"># 查找包含 TP53 的行</span>
grep -i <span class="hljs-string">"tp53"</span> genes.txt       <span class="hljs-comment"># 忽略大小写</span>
grep -v <span class="hljs-string">"comment"</span> file.txt     <span class="hljs-comment"># 反向匹配（排除）</span>
grep -c <span class="hljs-string">"gene"</span> file.txt        <span class="hljs-comment"># 计数</span>
grep -n <span class="hljs-string">"pattern"</span> file.txt     <span class="hljs-comment"># 显示行号</span>
grep -r <span class="hljs-string">"TODO"</span> src/            <span class="hljs-comment"># 递归搜索目录</span>

<span class="hljs-comment"># 管道组合</span>
ps aux | grep python           <span class="hljs-comment"># 找 python 进程</span>
<span class="hljs-built_in">cat</span> reads.fastq | grep -c <span class="hljs-string">"^@"</span> <span class="hljs-comment"># 统计 FASTQ 中序列条数</span>
</code></pre>
<h3>6.2 sed：流编辑</h3>
<pre><code class="hljs language-bash">sed -n <span class="hljs-string">'5,10p'</span> file.txt        <span class="hljs-comment"># 打印第 5-10 行</span>
sed <span class="hljs-string">'s/old/new/'</span> file.txt      <span class="hljs-comment"># 替换每行第一个匹配</span>
sed <span class="hljs-string">'s/old/new/g'</span> file.txt     <span class="hljs-comment"># 全局替换</span>
sed -i <span class="hljs-string">'s/old/new/g'</span> file.txt  <span class="hljs-comment"># 直接修改文件（-i 原地）</span>
sed <span class="hljs-string">'/^#/d'</span> config.conf        <span class="hljs-comment"># 删除注释行</span>
</code></pre>
<h3>6.3 awk：按列处理</h3>
<pre><code class="hljs language-bash">awk <span class="hljs-string">'{print $1}'</span> file.txt      <span class="hljs-comment"># 打印第一列</span>
awk -F<span class="hljs-string">','</span> <span class="hljs-string">'{print $1, $3}'</span> data.csv   <span class="hljs-comment"># 按逗号分隔</span>
awk <span class="hljs-string">'$3 > 50 {print $1}'</span> data.txt     <span class="hljs-comment"># 条件过滤</span>
awk <span class="hljs-string">'{sum += $2} END {print sum}'</span> data.txt   <span class="hljs-comment"># 第二列求和</span>
awk <span class="hljs-string">'NR > 1 {print}'</span> file.txt  <span class="hljs-comment"># 跳过表头</span>
</code></pre>
<blockquote>
<p>awk 字段分隔符默认空白；<code>-F','</code> 指定逗号。</p>
</blockquote>
<h2>7. 进程管理</h2>
<pre><code class="hljs language-bash">ps aux                <span class="hljs-comment"># 查看所有进程</span>
ps aux | grep python  <span class="hljs-comment"># 查找特定进程</span>
top                   <span class="hljs-comment"># 实时监控（q 退出）</span>
htop                  <span class="hljs-comment"># 增强版监控（更直观）</span>

<span class="hljs-built_in">kill</span> PID              <span class="hljs-comment"># 终止进程</span>
<span class="hljs-built_in">kill</span> -9 PID           <span class="hljs-comment"># 强制终止</span>
pkill python          <span class="hljs-comment"># 按名称终止</span>

<span class="hljs-comment"># 后台运行</span>
python script.py &#x26;    <span class="hljs-comment"># 后台运行</span>
<span class="hljs-built_in">nohup</span> python script.py > log.txt 2>&#x26;1 &#x26;   <span class="hljs-comment"># 脱离终端运行（服务器常用）</span>
<span class="hljs-built_in">jobs</span>                  <span class="hljs-comment"># 查看后台任务</span>
<span class="hljs-built_in">fg</span>                    <span class="hljs-comment"># 调回前台</span>
</code></pre>
<blockquote>
<p>远程服务器跑长任务用 <code>nohup ... &#x26;</code> 或 <code>tmux</code> / <code>screen</code>，避免终端断开导致任务中断。</p>
</blockquote>
<h2>8. 实用组合</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 统计文件数量与大小</span>
<span class="hljs-built_in">ls</span> -lh | awk <span class="hljs-string">'{print $5, $9}'</span>
<span class="hljs-built_in">du</span> -sh data/          <span class="hljs-comment"># 目录总大小</span>
<span class="hljs-built_in">df</span> -h                 <span class="hljs-comment"># 磁盘空间</span>

<span class="hljs-comment"># 压缩与解压</span>
tar -czf archive.tar.gz <span class="hljs-built_in">dir</span>/    <span class="hljs-comment"># 打包压缩</span>
tar -xzf archive.tar.gz         <span class="hljs-comment"># 解压</span>
zip -r archive.zip <span class="hljs-built_in">dir</span>/
unzip archive.zip

<span class="hljs-comment"># 查找文件</span>
find . -name <span class="hljs-string">"*.log"</span>            <span class="hljs-comment"># 按名查找</span>
find . -size +100M              <span class="hljs-comment"># 按大小查找</span>
<span class="hljs-built_in">which</span> python                    <span class="hljs-comment"># 命令所在路径</span>
</code></pre>
<h2>9. WSL 与开发环境</h2>
<p>Windows 下推荐 WSL（Windows Subsystem for Linux）：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 安装（管理员 PowerShell）</span>
wsl --install

<span class="hljs-comment"># 进入 Ubuntu 子系统</span>
wsl

<span class="hljs-comment"># 在 WSL 中配置开发环境</span>
<span class="hljs-built_in">sudo</span> apt update
<span class="hljs-built_in">sudo</span> apt install build-essential git python3 python3-pip
</code></pre>
<p>VS Code 连接 WSL：安装 Remote - WSL 扩展后，<code>code .</code> 直接打开。</p>
<h2>10. 小结</h2>
<ul>
<li>导航：<code>cd</code> / <code>ls</code> / <code>pwd</code>；操作：<code>cp</code> / <code>mv</code> / <code>rm</code> / <code>mkdir</code></li>
<li>三剑客：<code>grep</code> 搜索、<code>sed</code> 替换、<code>awk</code> 按列处理</li>
<li>管道 <code>|</code> 与重定向 <code>></code> 是组合命令的核心</li>
<li>权限 <code>chmod 755</code> / <code>chmod +x</code>；进程 <code>ps</code> / <code>kill</code> / <code>nohup</code></li>
</ul>
<p>下一篇将介绍 Bash 编程：把命令组织成可复用的脚本。</p>`,w_=`<h1>Python Advanced: Functions, Classes, and Modules</h1>
<p>After mastering the basics, this article dives into Python's functional and object-oriented programming, helping you write clean, reusable code.</p>
<h2>1. Functions</h2>
<h3>1.1 Definition and Call</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">greet</span>(<span class="hljs-params">name, greeting=<span class="hljs-string">"Hello"</span></span>):
    <span class="hljs-string">"""Return a greeting. greeting is a parameter with a default value."""</span>
    <span class="hljs-keyword">return</span> <span class="hljs-string">f"<span class="hljs-subst">{greeting}</span>, <span class="hljs-subst">{name}</span>!"</span>

<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>))            <span class="hljs-comment"># Hello, zorrooz!</span>
<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"Hi"</span>))      <span class="hljs-comment"># Hi, zorrooz!</span>
</code></pre>
<h3>1.2 Parameter Passing</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">show</span>(<span class="hljs-params">a, b, *args, **kwargs</span>):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"a=<span class="hljs-subst">{a}</span>, b=<span class="hljs-subst">{b}</span>"</span>)
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Positional args: <span class="hljs-subst">{args}</span>"</span>)      <span class="hljs-comment"># tuple</span>
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Keyword args: <span class="hljs-subst">{kwargs}</span>"</span>)  <span class="hljs-comment"># dict</span>

show(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, x=<span class="hljs-number">5</span>, y=<span class="hljs-number">6</span>)
<span class="hljs-comment"># a=1, b=2</span>
<span class="hljs-comment"># Positional args: (3, 4)</span>
<span class="hljs-comment"># Keyword args: {'x': 5, 'y': 6}</span>
</code></pre>
<ul>
<li><code>*args</code>: collects extra positional arguments into a tuple</li>
<li><code>**kwargs</code>: collects extra keyword arguments into a dictionary</li>
</ul>
<h3>1.3 Unpacking Arguments</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">add</span>(<span class="hljs-params">a, b, c</span>):
    <span class="hljs-keyword">return</span> a + b + c

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>]
<span class="hljs-built_in">print</span>(add(*nums))            <span class="hljs-comment"># 6, list unpacking</span>

data = {<span class="hljs-string">"a"</span>: <span class="hljs-number">1</span>, <span class="hljs-string">"b"</span>: <span class="hljs-number">2</span>, <span class="hljs-string">"c"</span>: <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(add(**data))           <span class="hljs-comment"># 6, dict unpacking</span>
</code></pre>
<h3>1.4 Lambda Anonymous Functions</h3>
<pre><code class="hljs language-python">square = <span class="hljs-keyword">lambda</span> x: x ** <span class="hljs-number">2</span>
<span class="hljs-built_in">print</span>(square(<span class="hljs-number">5</span>))             <span class="hljs-comment"># 25</span>

<span class="hljs-comment"># Combined with sorted / map / filter</span>
words = [<span class="hljs-string">"banana"</span>, <span class="hljs-string">"apple"</span>, <span class="hljs-string">"cherry"</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">sorted</span>(words, key=<span class="hljs-keyword">lambda</span> w: <span class="hljs-built_in">len</span>(w)))
<span class="hljs-comment"># ['apple', 'banana', 'cherry']</span>

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, <span class="hljs-number">5</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">map</span>(<span class="hljs-keyword">lambda</span> x: x * <span class="hljs-number">2</span>, nums)))       <span class="hljs-comment"># [2, 4, 6, 8, 10]</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">filter</span>(<span class="hljs-keyword">lambda</span> x: x % <span class="hljs-number">2</span> == <span class="hljs-number">0</span>, nums)))  <span class="hljs-comment"># [2, 4]</span>
</code></pre>
<h3>1.5 Closures</h3>
<p>Define a function inside another function and reference outer variables:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">make_counter</span>():
    count = <span class="hljs-number">0</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">counter</span>():
        <span class="hljs-keyword">nonlocal</span> count
        count += <span class="hljs-number">1</span>
        <span class="hljs-keyword">return</span> count
    <span class="hljs-keyword">return</span> counter

c = make_counter()
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 1</span>
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 2</span>
</code></pre>
<p><code>nonlocal</code> is used to modify variables in the outer function.</p>
<h2>2. Decorators</h2>
<p>A decorator is a higher-order function that "takes a function and returns a function," used to enhance its behavior without modifying the original function:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> time

<span class="hljs-keyword">def</span> <span class="hljs-title function_">timer</span>(<span class="hljs-params">func</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{func.__name__}</span> took <span class="hljs-subst">{elapsed:<span class="hljs-number">.4</span>f}</span>s"</span>)
        <span class="hljs-keyword">return</span> result
    <span class="hljs-keyword">return</span> wrapper

<span class="hljs-meta">@timer</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">slow_sum</span>(<span class="hljs-params">n</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">sum</span>(<span class="hljs-built_in">range</span>(n))

<span class="hljs-built_in">print</span>(slow_sum(<span class="hljs-number">10_000_000</span>))
<span class="hljs-comment"># slow_sum took 0.35xx s</span>
</code></pre>
<p>Decorator with parameters:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">repeat</span>(<span class="hljs-params">times</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">decorator</span>(<span class="hljs-params">func</span>):
        <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
            <span class="hljs-keyword">for</span> _ <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(times):
                func(*args, **kwargs)
        <span class="hljs-keyword">return</span> wrapper
    <span class="hljs-keyword">return</span> decorator

<span class="hljs-meta">@repeat(<span class="hljs-params"><span class="hljs-number">3</span></span>)</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">hello</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Hi!"</span>)

hello()
<span class="hljs-comment"># Hi! Hi! Hi!</span>
</code></pre>
<h2>3. Classes and Object-Oriented Programming</h2>
<h3>3.1 Basic Definition</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Sequence</span>:
    <span class="hljs-string">"""Represents a biological sequence."""</span>

    <span class="hljs-comment"># Class attribute: shared by all instances</span>
    alphabet = <span class="hljs-string">"ACGT"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-string">"""Constructor: initialize instance attributes."""</span>
        <span class="hljs-variable language_">self</span>.seq_id = seq_id
        <span class="hljs-variable language_">self</span>.seq = seq.upper()

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">length</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">gc_content</span>(<span class="hljs-params">self</span>):
        gc = <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"G"</span>) + <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"C"</span>)
        <span class="hljs-keyword">return</span> gc / <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq) * <span class="hljs-number">100</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Sequence(<span class="hljs-subst">{self.seq_id}</span>, <span class="hljs-subst">{self.length()}</span>bp)"</span>

s = <span class="hljs-type">Sequence</span>(<span class="hljs-string">"seq1"</span>, <span class="hljs-string">"atgcgta"</span>)
<span class="hljs-built_in">print</span>(s.length())        <span class="hljs-comment"># 7</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"GC: <span class="hljs-subst">{s.gc_content():<span class="hljs-number">.1</span>f}</span>%"</span>)   <span class="hljs-comment"># GC: 57.1%</span>
<span class="hljs-built_in">print</span>(s)                 <span class="hljs-comment"># Sequence(seq1, 7bp)</span>
</code></pre>
<h3>3.2 Inheritance</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Protein</span>(<span class="hljs-title class_ inherited__">Sequence</span>):
    alphabet = <span class="hljs-string">"ACDEFGHIKLMNPQRSTVWY"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-built_in">super</span>().__init__(seq_id, seq)   <span class="hljs-comment"># Call the parent class constructor</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">molecular_weight</span>(<span class="hljs-params">self</span>):
        <span class="hljs-comment"># Simplified calculation: each amino acid is approximately 110 Da</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.length() * <span class="hljs-number">110</span>

p = Protein(<span class="hljs-string">"prot1"</span>, <span class="hljs-string">"MKWVTFISLL"</span>)
<span class="hljs-built_in">print</span>(p.molecular_weight())   <span class="hljs-comment"># 1210</span>
</code></pre>
<h3>3.3 Magic Methods</h3>
<p>Common magic methods allow objects to support built-in operations:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Vector</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, x, y</span>):
        <span class="hljs-variable language_">self</span>.x = x
        <span class="hljs-variable language_">self</span>.y = y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__add__</span>(<span class="hljs-params">self, other</span>):      <span class="hljs-comment"># +</span>
        <span class="hljs-keyword">return</span> Vector(<span class="hljs-variable language_">self</span>.x + other.x, <span class="hljs-variable language_">self</span>.y + other.y)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__eq__</span>(<span class="hljs-params">self, other</span>):       <span class="hljs-comment"># ==</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.x == other.x <span class="hljs-keyword">and</span> <span class="hljs-variable language_">self</span>.y == other.y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):            <span class="hljs-comment"># print / repr</span>
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Vector(<span class="hljs-subst">{self.x}</span>, <span class="hljs-subst">{self.y}</span>)"</span>

v1 = Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>)
v2 = Vector(<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(v1 + v2)          <span class="hljs-comment"># Vector(4, 6)</span>
<span class="hljs-built_in">print</span>(v1 == Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># True</span>
</code></pre>
<h3>3.4 Properties (property)</h3>
<p>Use @property to turn methods into attribute access, and add validation:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Person</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, name</span>):
        <span class="hljs-variable language_">self</span>._name = name

<span class="hljs-meta">    @property</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>._name

<span class="hljs-meta">    @name.setter</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self, value</span>):
        <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> value.strip():
            <span class="hljs-keyword">raise</span> ValueError(<span class="hljs-string">"Name cannot be empty"</span>)
        <span class="hljs-variable language_">self</span>._name = value

p = Person(<span class="hljs-string">"zorrooz"</span>)
p.name = <span class="hljs-string">"bio"</span>
<span class="hljs-built_in">print</span>(p.name)
</code></pre>
<h2>4. Exception Handling</h2>
<pre><code class="hljs language-python"><span class="hljs-keyword">try</span>:
    num = <span class="hljs-built_in">int</span>(<span class="hljs-built_in">input</span>(<span class="hljs-string">"Enter an integer: "</span>))
    result = <span class="hljs-number">100</span> / num
<span class="hljs-keyword">except</span> ValueError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Not an integer"</span>)
<span class="hljs-keyword">except</span> ZeroDivisionError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Cannot divide by zero"</span>)
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Result: <span class="hljs-subst">{result}</span>"</span>)
<span class="hljs-keyword">finally</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"This always executes regardless of errors"</span>)
</code></pre>
<p>Custom exceptions:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">InvalidSequenceError</span>(<span class="hljs-title class_ inherited__">Exception</span>):
    <span class="hljs-keyword">pass</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">validate</span>(<span class="hljs-params">seq</span>):
    <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> <span class="hljs-built_in">set</span>(seq) &#x3C;= <span class="hljs-built_in">set</span>(<span class="hljs-string">"ACGT"</span>):
        <span class="hljs-keyword">raise</span> InvalidSequenceError(<span class="hljs-string">f"Contains invalid characters: <span class="hljs-subst">{seq}</span>"</span>)

<span class="hljs-keyword">try</span>:
    validate(<span class="hljs-string">"ATGXYZ"</span>)
<span class="hljs-keyword">except</span> InvalidSequenceError <span class="hljs-keyword">as</span> e:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"Validation failed: <span class="hljs-subst">{e}</span>"</span>)
</code></pre>
<h2>5. Modules and Packages</h2>
<h3>5.1 Module Imports</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> math                       <span class="hljs-comment"># Entire module</span>
<span class="hljs-keyword">from</span> math <span class="hljs-keyword">import</span> sqrt, pi         <span class="hljs-comment"># Import specific names</span>
<span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np                <span class="hljs-comment"># Alias</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter   <span class="hljs-comment"># Common: counting</span>

<span class="hljs-comment"># Counter example</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter
cnt = Counter(<span class="hljs-string">"ATGCCGA"</span>)
<span class="hljs-built_in">print</span>(cnt)          <span class="hljs-comment"># Counter({'C': 2, 'G': 2, 'A': 2, 'T': 1})</span>
<span class="hljs-built_in">print</span>(cnt.most_common(<span class="hljs-number">2</span>))
</code></pre>
<h3>5.2 Package Structure</h3>
<pre><code>myproject/
├── __init__.py        # Marks it as a package (optional in Python 3.3+)
├── utils/
│   ├── __init__.py
│   └── fasta.py       # Defines read_fasta()
└── main.py
</code></pre>
<pre><code class="hljs language-python"><span class="hljs-comment"># main.py</span>
<span class="hljs-keyword">from</span> utils.fasta <span class="hljs-keyword">import</span> read_fasta   <span class="hljs-comment"># Relative package import</span>
</code></pre>
<h3>5.3 <code>if __name__ == "__main__"</code></h3>
<p>Allows the script to be both imported and run directly:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">main</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Running main logic"</span>)

<span class="hljs-keyword">if</span> __name__ == <span class="hljs-string">"__main__"</span>:
    main()
</code></pre>
<h2>6. Type Hints (typing)</h2>
<p>Type hints improve readability and work with IDE static checks:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">from</span> typing <span class="hljs-keyword">import</span> <span class="hljs-type">List</span>, <span class="hljs-type">Dict</span>, <span class="hljs-type">Optional</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">count_bases</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Dict</span>[<span class="hljs-built_in">str</span>, <span class="hljs-built_in">int</span>]:
    <span class="hljs-string">"""Return the occurrence count of each base."""</span>
    <span class="hljs-keyword">return</span> {b: seq.count(b) <span class="hljs-keyword">for</span> b <span class="hljs-keyword">in</span> <span class="hljs-string">"ACGT"</span>}

<span class="hljs-keyword">def</span> <span class="hljs-title function_">find_motif</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span>, motif: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Optional</span>[<span class="hljs-built_in">int</span>]:
    idx = seq.find(motif)
    <span class="hljs-keyword">return</span> idx <span class="hljs-keyword">if</span> idx >= <span class="hljs-number">0</span> <span class="hljs-keyword">else</span> <span class="hljs-literal">None</span>

<span class="hljs-built_in">print</span>(count_bases(<span class="hljs-string">"ATGCCGA"</span>))
</code></pre>
<h2>7. Summary</h2>
<ul>
<li><code>*args</code> / <code>**kwargs</code>, lambda, closures, and decorators are core to functional programming</li>
<li>Classes: <code>__init__</code>, inheritance, magic methods, <code>@property</code></li>
<li>Exceptions: <code>try/except/else/finally</code>, custom exceptions inherit <code>Exception</code></li>
<li>Modularity: package directories + <code>if __name__ == "__main__"</code> guard</li>
<li>Type hints make code more maintainable</li>
</ul>
<p>The next article will cover Python data processing in practice: file I/O, regular expressions, and NumPy/Pandas.</p>`,C_=`<h1>Python 进阶：函数、类与模块</h1>
<p>掌握基础语法后，本文深入 Python 的函数式与面向对象编程，帮助你写出结构清晰、可复用的代码。</p>
<h2>1. 函数</h2>
<h3>1.1 定义与调用</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">greet</span>(<span class="hljs-params">name, greeting=<span class="hljs-string">"Hello"</span></span>):
    <span class="hljs-string">"""返回问候语。greeting 是带默认值的参数。"""</span>
    <span class="hljs-keyword">return</span> <span class="hljs-string">f"<span class="hljs-subst">{greeting}</span>, <span class="hljs-subst">{name}</span>!"</span>

<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>))            <span class="hljs-comment"># Hello, zorrooz!</span>
<span class="hljs-built_in">print</span>(greet(<span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"Hi"</span>))      <span class="hljs-comment"># Hi, zorrooz!</span>
</code></pre>
<h3>1.2 参数传递</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">show</span>(<span class="hljs-params">a, b, *args, **kwargs</span>):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"a=<span class="hljs-subst">{a}</span>, b=<span class="hljs-subst">{b}</span>"</span>)
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"位置参数: <span class="hljs-subst">{args}</span>"</span>)      <span class="hljs-comment"># 元组</span>
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"关键字参数: <span class="hljs-subst">{kwargs}</span>"</span>)  <span class="hljs-comment"># 字典</span>

show(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, x=<span class="hljs-number">5</span>, y=<span class="hljs-number">6</span>)
<span class="hljs-comment"># a=1, b=2</span>
<span class="hljs-comment"># 位置参数: (3, 4)</span>
<span class="hljs-comment"># 关键字参数: {'x': 5, 'y': 6}</span>
</code></pre>
<ul>
<li><code>*args</code>：收集多余的位置参数为元组</li>
<li><code>**kwargs</code>：收集多余的关键字参数为字典</li>
</ul>
<h3>1.3 解包传参</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">add</span>(<span class="hljs-params">a, b, c</span>):
    <span class="hljs-keyword">return</span> a + b + c

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>]
<span class="hljs-built_in">print</span>(add(*nums))            <span class="hljs-comment"># 6，列表解包</span>

data = {<span class="hljs-string">"a"</span>: <span class="hljs-number">1</span>, <span class="hljs-string">"b"</span>: <span class="hljs-number">2</span>, <span class="hljs-string">"c"</span>: <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(add(**data))           <span class="hljs-comment"># 6，字典解包</span>
</code></pre>
<h3>1.4 lambda 匿名函数</h3>
<pre><code class="hljs language-python">square = <span class="hljs-keyword">lambda</span> x: x ** <span class="hljs-number">2</span>
<span class="hljs-built_in">print</span>(square(<span class="hljs-number">5</span>))             <span class="hljs-comment"># 25</span>

<span class="hljs-comment"># 与 sorted / map / filter 配合</span>
words = [<span class="hljs-string">"banana"</span>, <span class="hljs-string">"apple"</span>, <span class="hljs-string">"cherry"</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">sorted</span>(words, key=<span class="hljs-keyword">lambda</span> w: <span class="hljs-built_in">len</span>(w)))
<span class="hljs-comment"># ['apple', 'banana', 'cherry']</span>

nums = [<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>, <span class="hljs-number">5</span>]
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">map</span>(<span class="hljs-keyword">lambda</span> x: x * <span class="hljs-number">2</span>, nums)))       <span class="hljs-comment"># [2, 4, 6, 8, 10]</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">list</span>(<span class="hljs-built_in">filter</span>(<span class="hljs-keyword">lambda</span> x: x % <span class="hljs-number">2</span> == <span class="hljs-number">0</span>, nums)))  <span class="hljs-comment"># [2, 4]</span>
</code></pre>
<h3>1.5 闭包</h3>
<p>函数内部定义函数并引用外部变量：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">make_counter</span>():
    count = <span class="hljs-number">0</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">counter</span>():
        <span class="hljs-keyword">nonlocal</span> count
        count += <span class="hljs-number">1</span>
        <span class="hljs-keyword">return</span> count
    <span class="hljs-keyword">return</span> counter

c = make_counter()
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 1</span>
<span class="hljs-built_in">print</span>(c())   <span class="hljs-comment"># 2</span>
</code></pre>
<p><code>nonlocal</code> 用于修改外层函数中的变量。</p>
<h2>2. 装饰器</h2>
<p>装饰器是"接收函数、返回函数"的高阶函数，用于在不修改原函数的情况下增强其行为：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> time

<span class="hljs-keyword">def</span> <span class="hljs-title function_">timer</span>(<span class="hljs-params">func</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
        start = time.perf_counter()
        result = func(*args, **kwargs)
        elapsed = time.perf_counter() - start
        <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{func.__name__}</span> 耗时 <span class="hljs-subst">{elapsed:<span class="hljs-number">.4</span>f}</span>s"</span>)
        <span class="hljs-keyword">return</span> result
    <span class="hljs-keyword">return</span> wrapper

<span class="hljs-meta">@timer</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">slow_sum</span>(<span class="hljs-params">n</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">sum</span>(<span class="hljs-built_in">range</span>(n))

<span class="hljs-built_in">print</span>(slow_sum(<span class="hljs-number">10_000_000</span>))
<span class="hljs-comment"># slow_sum 耗时 0.35xx s</span>
</code></pre>
<p>带参数的装饰器：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">repeat</span>(<span class="hljs-params">times</span>):
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">decorator</span>(<span class="hljs-params">func</span>):
        <span class="hljs-keyword">def</span> <span class="hljs-title function_">wrapper</span>(<span class="hljs-params">*args, **kwargs</span>):
            <span class="hljs-keyword">for</span> _ <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(times):
                func(*args, **kwargs)
        <span class="hljs-keyword">return</span> wrapper
    <span class="hljs-keyword">return</span> decorator

<span class="hljs-meta">@repeat(<span class="hljs-params"><span class="hljs-number">3</span></span>)</span>
<span class="hljs-keyword">def</span> <span class="hljs-title function_">hello</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"Hi!"</span>)

hello()
<span class="hljs-comment"># Hi! Hi! Hi!</span>
</code></pre>
<h2>3. 类与面向对象</h2>
<h3>3.1 基本定义</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Sequence</span>:
    <span class="hljs-string">"""表示一条生物序列。"""</span>

    <span class="hljs-comment"># 类属性：所有实例共享</span>
    alphabet = <span class="hljs-string">"ACGT"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-string">"""构造方法：初始化实例属性。"""</span>
        <span class="hljs-variable language_">self</span>.seq_id = seq_id
        <span class="hljs-variable language_">self</span>.seq = seq.upper()

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">length</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">gc_content</span>(<span class="hljs-params">self</span>):
        gc = <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"G"</span>) + <span class="hljs-variable language_">self</span>.seq.count(<span class="hljs-string">"C"</span>)
        <span class="hljs-keyword">return</span> gc / <span class="hljs-built_in">len</span>(<span class="hljs-variable language_">self</span>.seq) * <span class="hljs-number">100</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Sequence(<span class="hljs-subst">{self.seq_id}</span>, <span class="hljs-subst">{self.length()}</span>bp)"</span>

s = <span class="hljs-type">Sequence</span>(<span class="hljs-string">"seq1"</span>, <span class="hljs-string">"atgcgta"</span>)
<span class="hljs-built_in">print</span>(s.length())        <span class="hljs-comment"># 7</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"GC: <span class="hljs-subst">{s.gc_content():<span class="hljs-number">.1</span>f}</span>%"</span>)   <span class="hljs-comment"># GC: 57.1%</span>
<span class="hljs-built_in">print</span>(s)                 <span class="hljs-comment"># Sequence(seq1, 7bp)</span>
</code></pre>
<h3>3.2 继承</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Protein</span>(<span class="hljs-title class_ inherited__">Sequence</span>):
    alphabet = <span class="hljs-string">"ACDEFGHIKLMNPQRSTVWY"</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, seq_id, seq</span>):
        <span class="hljs-built_in">super</span>().__init__(seq_id, seq)   <span class="hljs-comment"># 调用父类构造方法</span>

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">molecular_weight</span>(<span class="hljs-params">self</span>):
        <span class="hljs-comment"># 简化计算：每个氨基酸约 110 Da</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.length() * <span class="hljs-number">110</span>

p = Protein(<span class="hljs-string">"prot1"</span>, <span class="hljs-string">"MKWVTFISLL"</span>)
<span class="hljs-built_in">print</span>(p.molecular_weight())   <span class="hljs-comment"># 1210</span>
</code></pre>
<h3>3.3 魔术方法</h3>
<p>常用魔术方法让对象支持内置操作：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Vector</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, x, y</span>):
        <span class="hljs-variable language_">self</span>.x = x
        <span class="hljs-variable language_">self</span>.y = y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__add__</span>(<span class="hljs-params">self, other</span>):      <span class="hljs-comment"># +</span>
        <span class="hljs-keyword">return</span> Vector(<span class="hljs-variable language_">self</span>.x + other.x, <span class="hljs-variable language_">self</span>.y + other.y)

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__eq__</span>(<span class="hljs-params">self, other</span>):       <span class="hljs-comment"># ==</span>
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>.x == other.x <span class="hljs-keyword">and</span> <span class="hljs-variable language_">self</span>.y == other.y

    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__repr__</span>(<span class="hljs-params">self</span>):            <span class="hljs-comment"># print / repr</span>
        <span class="hljs-keyword">return</span> <span class="hljs-string">f"Vector(<span class="hljs-subst">{self.x}</span>, <span class="hljs-subst">{self.y}</span>)"</span>

v1 = Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>)
v2 = Vector(<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(v1 + v2)          <span class="hljs-comment"># Vector(4, 6)</span>
<span class="hljs-built_in">print</span>(v1 == Vector(<span class="hljs-number">1</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># True</span>
</code></pre>
<h3>3.4 属性（property）</h3>
<p>用 <code>@property</code> 把方法变成属性访问，可加入校验：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">Person</span>:
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">__init__</span>(<span class="hljs-params">self, name</span>):
        <span class="hljs-variable language_">self</span>._name = name

<span class="hljs-meta">    @property</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self</span>):
        <span class="hljs-keyword">return</span> <span class="hljs-variable language_">self</span>._name

<span class="hljs-meta">    @name.setter</span>
    <span class="hljs-keyword">def</span> <span class="hljs-title function_">name</span>(<span class="hljs-params">self, value</span>):
        <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> value.strip():
            <span class="hljs-keyword">raise</span> ValueError(<span class="hljs-string">"名字不能为空"</span>)
        <span class="hljs-variable language_">self</span>._name = value

p = Person(<span class="hljs-string">"zorrooz"</span>)
p.name = <span class="hljs-string">"bio"</span>
<span class="hljs-built_in">print</span>(p.name)
</code></pre>
<h2>4. 异常处理</h2>
<pre><code class="hljs language-python"><span class="hljs-keyword">try</span>:
    num = <span class="hljs-built_in">int</span>(<span class="hljs-built_in">input</span>(<span class="hljs-string">"输入一个整数："</span>))
    result = <span class="hljs-number">100</span> / num
<span class="hljs-keyword">except</span> ValueError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"输入的不是整数"</span>)
<span class="hljs-keyword">except</span> ZeroDivisionError:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"不能除以零"</span>)
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"结果: <span class="hljs-subst">{result}</span>"</span>)
<span class="hljs-keyword">finally</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"无论是否出错都会执行"</span>)
</code></pre>
<p>自定义异常：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">class</span> <span class="hljs-title class_">InvalidSequenceError</span>(<span class="hljs-title class_ inherited__">Exception</span>):
    <span class="hljs-keyword">pass</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">validate</span>(<span class="hljs-params">seq</span>):
    <span class="hljs-keyword">if</span> <span class="hljs-keyword">not</span> <span class="hljs-built_in">set</span>(seq) &#x3C;= <span class="hljs-built_in">set</span>(<span class="hljs-string">"ACGT"</span>):
        <span class="hljs-keyword">raise</span> InvalidSequenceError(<span class="hljs-string">f"包含非法字符: <span class="hljs-subst">{seq}</span>"</span>)

<span class="hljs-keyword">try</span>:
    validate(<span class="hljs-string">"ATGXYZ"</span>)
<span class="hljs-keyword">except</span> InvalidSequenceError <span class="hljs-keyword">as</span> e:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"校验失败: <span class="hljs-subst">{e}</span>"</span>)
</code></pre>
<h2>5. 模块与包</h2>
<h3>5.1 模块导入</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> math                       <span class="hljs-comment"># 整个模块</span>
<span class="hljs-keyword">from</span> math <span class="hljs-keyword">import</span> sqrt, pi         <span class="hljs-comment"># 导入指定名字</span>
<span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np                <span class="hljs-comment"># 别名</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter   <span class="hljs-comment"># 常用：计数</span>

<span class="hljs-comment"># Counter 示例</span>
<span class="hljs-keyword">from</span> collections <span class="hljs-keyword">import</span> Counter
cnt = Counter(<span class="hljs-string">"ATGCCGA"</span>)
<span class="hljs-built_in">print</span>(cnt)          <span class="hljs-comment"># Counter({'C': 2, 'G': 2, 'A': 2, 'T': 1})</span>
<span class="hljs-built_in">print</span>(cnt.most_common(<span class="hljs-number">2</span>))
</code></pre>
<h3>5.2 包结构</h3>
<pre><code>myproject/
├── __init__.py        # 标记为包（Python 3.3+ 可省略）
├── utils/
│   ├── __init__.py
│   └── fasta.py       # 定义 read_fasta()
└── main.py
</code></pre>
<pre><code class="hljs language-python"><span class="hljs-comment"># main.py</span>
<span class="hljs-keyword">from</span> utils.fasta <span class="hljs-keyword">import</span> read_fasta   <span class="hljs-comment"># 相对包路径导入</span>
</code></pre>
<h3>5.3 <code>if __name__ == "__main__"</code></h3>
<p>让脚本既可被导入也可直接运行：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">main</span>():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"运行主逻辑"</span>)

<span class="hljs-keyword">if</span> __name__ == <span class="hljs-string">"__main__"</span>:
    main()
</code></pre>
<h2>6. 类型提示（typing）</h2>
<p>类型提示提高可读性，配合 IDE 静态检查：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">from</span> typing <span class="hljs-keyword">import</span> <span class="hljs-type">List</span>, <span class="hljs-type">Dict</span>, <span class="hljs-type">Optional</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">count_bases</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Dict</span>[<span class="hljs-built_in">str</span>, <span class="hljs-built_in">int</span>]:
    <span class="hljs-string">"""返回每种碱基的出现次数。"""</span>
    <span class="hljs-keyword">return</span> {b: seq.count(b) <span class="hljs-keyword">for</span> b <span class="hljs-keyword">in</span> <span class="hljs-string">"ACGT"</span>}

<span class="hljs-keyword">def</span> <span class="hljs-title function_">find_motif</span>(<span class="hljs-params">seq: <span class="hljs-built_in">str</span>, motif: <span class="hljs-built_in">str</span></span>) -> <span class="hljs-type">Optional</span>[<span class="hljs-built_in">int</span>]:
    idx = seq.find(motif)
    <span class="hljs-keyword">return</span> idx <span class="hljs-keyword">if</span> idx >= <span class="hljs-number">0</span> <span class="hljs-keyword">else</span> <span class="hljs-literal">None</span>

<span class="hljs-built_in">print</span>(count_bases(<span class="hljs-string">"ATGCCGA"</span>))
</code></pre>
<h2>7. 小结</h2>
<ul>
<li><code>*args</code> / <code>**kwargs</code>、lambda、闭包、装饰器是函数式编程核心</li>
<li>类：<code>__init__</code>、继承、魔术方法、<code>@property</code></li>
<li>异常：<code>try/except/else/finally</code>，自定义异常继承 <code>Exception</code></li>
<li>模块化：包目录 + <code>if __name__ == "__main__"</code> 守卫</li>
<li>类型提示让代码更可维护</li>
</ul>
<p>下一篇将介绍 Python 数据处理实战：文件 IO、正则表达式与 NumPy/Pandas。</p>`,k_=`<h1>Introduction to Python Programming: Environment, Syntax, and Data Types</h1>
<p>Python is an interpreted language with clean syntax and a rich ecosystem, widely used in scientific computing, data processing, and script automation. This article starts from scratch, covering environment setup, basic syntax, built-in data types, and control flow. All examples are directly runnable.</p>
<h2>1. Environment Setup</h2>
<h3>1.1 Installing Python</h3>
<p>Visit <a href="https://www.python.org/">python.org</a> to download the installer for your platform. When installing, be sure to check <strong>Add Python to PATH</strong>.</p>
<p>Verify the installation:</p>
<pre><code class="hljs language-bash">python --version
<span class="hljs-comment"># Python 3.12.x</span>
</code></pre>
<blockquote>
<p>On Linux / WSL, the <code>python3</code> command is typically used; on Windows, it is <code>python</code>.</p>
</blockquote>
<h3>1.2 Virtual Environment (Highly Recommended)</h3>
<p>A virtual environment isolates dependencies for each project, avoiding pollution of the global environment:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Create a virtual environment</span>
python -m venv .venv

<span class="hljs-comment"># Activate (Windows)</span>
.venv\\Scripts\\activate
<span class="hljs-comment"># Activate (Linux / macOS / WSL)</span>
<span class="hljs-built_in">source</span> .venv/bin/activate
</code></pre>
<p>After activation, <code>(.venv)</code> appears before the prompt. Install dependencies:</p>
<pre><code class="hljs language-bash">pip install numpy pandas
</code></pre>
<p>Export the dependency list to a file for easy reproduction:</p>
<pre><code class="hljs language-bash">pip freeze > requirements.txt
</code></pre>
<h3>1.3 Interactive Interpreter and Scripts</h3>
<pre><code class="hljs language-bash">python            <span class="hljs-comment"># Enter REPL interactive mode</span>
python script.py  <span class="hljs-comment"># Run a script file</span>
</code></pre>
<p>In the REPL, you can directly enter expressions and get immediate results, which is great for quickly validating ideas.</p>
<h2>2. First Program</h2>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(<span class="hljs-string">"Hello, Python!"</span>)
</code></pre>
<p>Save it as <code>hello.py</code> and run:</p>
<pre><code class="hljs language-bash">python hello.py
<span class="hljs-comment"># Hello, Python!</span>
</code></pre>
<h2>3. Basic Syntax</h2>
<h3>3.1 Comments</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Single-line comment</span>

<span class="hljs-string">"""
Multi-line comment (actually a multi-line string,
often used as documentation)
"""</span>
</code></pre>
<h3>3.2 Variables and Assignment</h3>
<p>Python is a dynamically typed language; variables do not require type declarations:</p>
<pre><code class="hljs language-python">name = <span class="hljs-string">"zorrooz"</span>     <span class="hljs-comment"># string</span>
age = <span class="hljs-number">25</span>             <span class="hljs-comment"># integer</span>
height = <span class="hljs-number">1.78</span>        <span class="hljs-comment"># float</span>
is_student = <span class="hljs-literal">True</span>    <span class="hljs-comment"># boolean</span>

<span class="hljs-comment"># Assign multiple variables at once</span>
a, b, c = <span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>

<span class="hljs-comment"># Swap two variables</span>
a, b = b, a
</code></pre>
<p>Variable naming follows PEP 8: lowercase letters + underscores (<code>snake_case</code>).</p>
<h3>3.3 Input and Output</h3>
<pre><code class="hljs language-python">name = <span class="hljs-built_in">input</span>(<span class="hljs-string">"请输入你的名字："</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，"</span>, name)

<span class="hljs-comment"># f-string formatting (Python 3.6+, recommended)</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"你好，<span class="hljs-subst">{name}</span>，今年 <span class="hljs-subst">{age}</span> 岁"</span>)

<span class="hljs-comment"># format method</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，{}，今年 {} 岁"</span>.<span class="hljs-built_in">format</span>(name, age))
</code></pre>
<h2>4. Built-in Data Types</h2>
<h3>4.1 Numbers</h3>
<pre><code class="hljs language-python">x = <span class="hljs-number">10</span>          <span class="hljs-comment"># int</span>
y = <span class="hljs-number">3.14</span>        <span class="hljs-comment"># float</span>
z = <span class="hljs-number">2</span> + <span class="hljs-number">3j</span>      <span class="hljs-comment"># complex</span>

<span class="hljs-comment"># Common operations</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> // <span class="hljs-number">2</span>)   <span class="hljs-comment"># integer division, result 3</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> % <span class="hljs-number">2</span>)    <span class="hljs-comment"># modulo, result 1</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">2</span> ** <span class="hljs-number">10</span>)  <span class="hljs-comment"># exponentiation, result 1024</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">round</span>(<span class="hljs-number">3.14159</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># round to 2 decimal places</span>
</code></pre>
<h3>4.2 Strings</h3>
<pre><code class="hljs language-python">s = <span class="hljs-string">"hello"</span>
t = <span class="hljs-string">'world'</span>

<span class="hljs-comment"># Concatenation and repetition</span>
<span class="hljs-built_in">print</span>(s + <span class="hljs-string">" "</span> + t)   <span class="hljs-comment"># hello world</span>
<span class="hljs-built_in">print</span>(s * <span class="hljs-number">3</span>)         <span class="hljs-comment"># hellohellohello</span>

<span class="hljs-comment"># Indexing and slicing (left-inclusive, right-exclusive)</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">0</span>])          <span class="hljs-comment"># h</span>
<span class="hljs-built_in">print</span>(s[-<span class="hljs-number">1</span>])         <span class="hljs-comment"># o</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])        <span class="hljs-comment"># el</span>

<span class="hljs-comment"># Common methods</span>
<span class="hljs-built_in">print</span>(s.upper())          <span class="hljs-comment"># HELLO</span>
<span class="hljs-built_in">print</span>(s.replace(<span class="hljs-string">"l"</span>, <span class="hljs-string">"L"</span>))  <span class="hljs-comment"># heLLo</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">len</span>(s))             <span class="hljs-comment"># 5</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">", "</span>.join([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>]))  <span class="hljs-comment"># a, b, c</span>
</code></pre>
<h3>4.3 Lists (list)</h3>
<p>Lists are ordered, mutable containers:</p>
<pre><code class="hljs language-python">fruits = [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>, <span class="hljs-string">"cherry"</span>]

fruits.append(<span class="hljs-string">"orange"</span>)      <span class="hljs-comment"># append at the end</span>
fruits.insert(<span class="hljs-number">0</span>, <span class="hljs-string">"grape"</span>)    <span class="hljs-comment"># insert at a specific position</span>
fruits.remove(<span class="hljs-string">"apple"</span>)       <span class="hljs-comment"># remove by value</span>
popped = fruits.pop()        <span class="hljs-comment"># pop the last element</span>

<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">0</span>])             <span class="hljs-comment"># indexing</span>
<span class="hljs-built_in">print</span>(fruits[-<span class="hljs-number">1</span>])            <span class="hljs-comment"># last element</span>
<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])           <span class="hljs-comment"># slicing</span>

<span class="hljs-comment"># List comprehension (very commonly used)</span>
squares = [x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>)]
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]</span>
</code></pre>
<h3>4.4 Tuples (tuple)</h3>
<p>Tuples are immutable lists, often used for returning multiple values:</p>
<pre><code class="hljs language-python">point = (<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
x, y = point          <span class="hljs-comment"># unpacking</span>
<span class="hljs-built_in">print</span>(x, y)           <span class="hljs-comment"># 3 4</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">minmax</span>(<span class="hljs-params">nums</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">min</span>(nums), <span class="hljs-built_in">max</span>(nums)

lo, hi = minmax([<span class="hljs-number">3</span>, <span class="hljs-number">1</span>, <span class="hljs-number">4</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>])
<span class="hljs-built_in">print</span>(lo, hi)         <span class="hljs-comment"># 1 5</span>
</code></pre>
<h3>4.5 Dictionaries (dict)</h3>
<p>Dictionaries are key-value containers with efficient lookup:</p>
<pre><code class="hljs language-python">person = {<span class="hljs-string">"name"</span>: <span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"age"</span>: <span class="hljs-number">25</span>}

<span class="hljs-built_in">print</span>(person[<span class="hljs-string">"name"</span>])            <span class="hljs-comment"># access value</span>
<span class="hljs-built_in">print</span>(person.get(<span class="hljs-string">"city"</span>, <span class="hljs-string">"未知"</span>))  <span class="hljs-comment"># safe access with default value</span>
person[<span class="hljs-string">"city"</span>] = <span class="hljs-string">"Lanzhou"</span>       <span class="hljs-comment"># add/update</span>

<span class="hljs-keyword">for</span> key, value <span class="hljs-keyword">in</span> person.items():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{key}</span>: <span class="hljs-subst">{value}</span>"</span>)

<span class="hljs-comment"># Dictionary comprehension</span>
squares = {x: x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>)}
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}</span>
</code></pre>
<h3>4.6 Sets (set)</h3>
<p>Sets are unordered, deduplicating containers:</p>
<pre><code class="hljs language-python">a = {<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(a)              <span class="hljs-comment"># {1, 2, 3}, automatically deduplicated</span>

b = {<span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>}
<span class="hljs-built_in">print</span>(a &#x26; b)          <span class="hljs-comment"># intersection {2, 3}</span>
<span class="hljs-built_in">print</span>(a | b)          <span class="hljs-comment"># union {1, 2, 3, 4}</span>
<span class="hljs-built_in">print</span>(a - b)          <span class="hljs-comment"># difference {1}</span>
</code></pre>
<h2>5. Control Flow</h2>
<h3>5.1 Conditional Statements</h3>
<pre><code class="hljs language-python">score = <span class="hljs-number">85</span>

<span class="hljs-keyword">if</span> score >= <span class="hljs-number">90</span>:
    grade = <span class="hljs-string">"A"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">80</span>:
    grade = <span class="hljs-string">"B"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">70</span>:
    grade = <span class="hljs-string">"C"</span>
<span class="hljs-keyword">else</span>:
    grade = <span class="hljs-string">"D"</span>

<span class="hljs-built_in">print</span>(<span class="hljs-string">f"等级：<span class="hljs-subst">{grade}</span>"</span>)

<span class="hljs-comment"># Ternary expression</span>
status = <span class="hljs-string">"pass"</span> <span class="hljs-keyword">if</span> score >= <span class="hljs-number">60</span> <span class="hljs-keyword">else</span> <span class="hljs-string">"fail"</span>
</code></pre>
<h3>5.2 for Loops</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Iterate over iterable objects</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>]:
    <span class="hljs-built_in">print</span>(fruit)

<span class="hljs-comment"># Generate number sequences with range</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):        <span class="hljs-comment"># 0 to 4</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">1</span>, <span class="hljs-number">10</span>, <span class="hljs-number">2</span>):  <span class="hljs-comment"># 1, 3, 5, 7, 9</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-comment"># Use enumerate to get index and value simultaneously</span>
<span class="hljs-keyword">for</span> idx, fruit <span class="hljs-keyword">in</span> <span class="hljs-built_in">enumerate</span>([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>]):
    <span class="hljs-built_in">print</span>(idx, fruit)
</code></pre>
<h3>5.3 while Loops</h3>
<pre><code class="hljs language-python">n = <span class="hljs-number">0</span>
<span class="hljs-keyword">while</span> n &#x3C; <span class="hljs-number">5</span>:
    <span class="hljs-built_in">print</span>(n)
    n += <span class="hljs-number">1</span>          <span class="hljs-comment"># Note: forgetting to increment causes an infinite loop</span>
</code></pre>
<h3>5.4 break / continue / else</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># break: terminate the loop early</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">3</span>:
        <span class="hljs-keyword">break</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># outputs 0 1 2</span>

<span class="hljs-comment"># continue: skip the current iteration</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">2</span>:
        <span class="hljs-keyword">continue</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># outputs 0 1 3 4</span>

<span class="hljs-comment"># for-else: executes when the loop completes normally (not broken)</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">3</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">99</span>:
        <span class="hljs-keyword">break</span>
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"循环正常结束"</span>)
</code></pre>
<h2>6. Comprehensive Exercise</h2>
<p>Count the occurrences of each word in a text:</p>
<pre><code class="hljs language-python">text = <span class="hljs-string">"python is fun and python is powerful"</span>

words = text.split()
counter = {}

<span class="hljs-keyword">for</span> word <span class="hljs-keyword">in</span> words:
    counter[word] = counter.get(word, <span class="hljs-number">0</span>) + <span class="hljs-number">1</span>

<span class="hljs-keyword">for</span> word, count <span class="hljs-keyword">in</span> <span class="hljs-built_in">sorted</span>(counter.items(), key=<span class="hljs-keyword">lambda</span> x: -x[<span class="hljs-number">1</span>]):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{word}</span>: <span class="hljs-subst">{count}</span>"</span>)
</code></pre>
<p>Output:</p>
<pre><code>python: 2
is: 2
fun: 1
and: 1
powerful: 1
</code></pre>
<h2>7. Summary</h2>
<ul>
<li>Use <code>venv</code> to isolate project dependencies, use <code>pip freeze</code> to record dependencies</li>
<li>Built-in types: <code>int</code> / <code>float</code> / <code>str</code> / <code>list</code> / <code>tuple</code> / <code>dict</code> / <code>set</code></li>
<li>List and dictionary comprehensions and f-strings are high-frequency syntax in daily use</li>
<li>Control flow: <code>if/elif/else</code>, <code>for</code>, <code>while</code>, <code>break</code> / <code>continue</code></li>
</ul>
<p>The next article will introduce functions, classes, and modules, moving into real engineering-style programming.</p>`,x_=`<h1>Python 编程入门：环境、语法与数据类型</h1>
<p>Python 是一门语法简洁、生态丰富的解释型语言，广泛应用于科学计算、数据处理与脚本自动化。本文从零开始，覆盖环境搭建、基础语法、内置数据类型与流程控制，所有示例均可直接运行。</p>
<h2>1. 环境搭建</h2>
<h3>1.1 安装 Python</h3>
<p>访问 <a href="https://www.python.org/">python.org</a> 下载对应平台的安装包。安装时务必勾选 <strong>Add Python to PATH</strong>。</p>
<p>验证安装：</p>
<pre><code class="hljs language-bash">python --version
<span class="hljs-comment"># Python 3.12.x</span>
</code></pre>
<blockquote>
<p>在 Linux / WSL 下通常使用 <code>python3</code> 命令，Windows 下为 <code>python</code>。</p>
</blockquote>
<h3>1.2 虚拟环境（强烈推荐）</h3>
<p>虚拟环境为每个项目隔离依赖，避免污染全局环境：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 创建虚拟环境</span>
python -m venv .venv

<span class="hljs-comment"># 激活（Windows）</span>
.venv\\Scripts\\activate
<span class="hljs-comment"># 激活（Linux / macOS / WSL）</span>
<span class="hljs-built_in">source</span> .venv/bin/activate
</code></pre>
<p>激活后提示符前会出现 <code>(.venv)</code>。安装依赖：</p>
<pre><code class="hljs language-bash">pip install numpy pandas
</code></pre>
<p>将依赖清单导出到文件，方便复现：</p>
<pre><code class="hljs language-bash">pip freeze > requirements.txt
</code></pre>
<h3>1.3 交互式解释器与脚本</h3>
<pre><code class="hljs language-bash">python            <span class="hljs-comment"># 进入 REPL 交互模式</span>
python script.py  <span class="hljs-comment"># 运行脚本文件</span>
</code></pre>
<p>在 REPL 中可直接输入表达式并立即得到结果，适合快速验证想法。</p>
<h2>2. 第一个程序</h2>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(<span class="hljs-string">"Hello, Python!"</span>)
</code></pre>
<p>保存为 <code>hello.py</code> 后运行：</p>
<pre><code class="hljs language-bash">python hello.py
<span class="hljs-comment"># Hello, Python!</span>
</code></pre>
<h2>3. 基础语法</h2>
<h3>3.1 注释</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 单行注释</span>

<span class="hljs-string">"""
多行注释（实际上是多行字符串，
但常被用作文档说明）
"""</span>
</code></pre>
<h3>3.2 变量与赋值</h3>
<p>Python 是动态类型语言，变量无需声明类型：</p>
<pre><code class="hljs language-python">name = <span class="hljs-string">"zorrooz"</span>     <span class="hljs-comment"># 字符串</span>
age = <span class="hljs-number">25</span>             <span class="hljs-comment"># 整数</span>
height = <span class="hljs-number">1.78</span>        <span class="hljs-comment"># 浮点数</span>
is_student = <span class="hljs-literal">True</span>    <span class="hljs-comment"># 布尔值</span>

<span class="hljs-comment"># 同时赋值多个变量</span>
a, b, c = <span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>

<span class="hljs-comment"># 交换两个变量的值</span>
a, b = b, a
</code></pre>
<p>变量命名遵循 PEP 8：小写字母 + 下划线（<code>snake_case</code>）。</p>
<h3>3.3 输入输出</h3>
<pre><code class="hljs language-python">name = <span class="hljs-built_in">input</span>(<span class="hljs-string">"请输入你的名字："</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，"</span>, name)

<span class="hljs-comment"># f-string 格式化（Python 3.6+，最推荐）</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"你好，<span class="hljs-subst">{name}</span>，今年 <span class="hljs-subst">{age}</span> 岁"</span>)

<span class="hljs-comment"># format 方法</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">"你好，{}，今年 {} 岁"</span>.<span class="hljs-built_in">format</span>(name, age))
</code></pre>
<h2>4. 内置数据类型</h2>
<h3>4.1 数字</h3>
<pre><code class="hljs language-python">x = <span class="hljs-number">10</span>          <span class="hljs-comment"># int</span>
y = <span class="hljs-number">3.14</span>        <span class="hljs-comment"># float</span>
z = <span class="hljs-number">2</span> + <span class="hljs-number">3j</span>      <span class="hljs-comment"># complex</span>

<span class="hljs-comment"># 常用运算</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> // <span class="hljs-number">2</span>)   <span class="hljs-comment"># 整除，得 3</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">7</span> % <span class="hljs-number">2</span>)    <span class="hljs-comment"># 取余，得 1</span>
<span class="hljs-built_in">print</span>(<span class="hljs-number">2</span> ** <span class="hljs-number">10</span>)  <span class="hljs-comment"># 幂运算，得 1024</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">round</span>(<span class="hljs-number">3.14159</span>, <span class="hljs-number">2</span>))  <span class="hljs-comment"># 四舍五入保留 2 位</span>
</code></pre>
<h3>4.2 字符串</h3>
<pre><code class="hljs language-python">s = <span class="hljs-string">"hello"</span>
t = <span class="hljs-string">'world'</span>

<span class="hljs-comment"># 拼接与重复</span>
<span class="hljs-built_in">print</span>(s + <span class="hljs-string">" "</span> + t)   <span class="hljs-comment"># hello world</span>
<span class="hljs-built_in">print</span>(s * <span class="hljs-number">3</span>)         <span class="hljs-comment"># hellohellohello</span>

<span class="hljs-comment"># 索引与切片（左闭右开）</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">0</span>])          <span class="hljs-comment"># h</span>
<span class="hljs-built_in">print</span>(s[-<span class="hljs-number">1</span>])         <span class="hljs-comment"># o</span>
<span class="hljs-built_in">print</span>(s[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])        <span class="hljs-comment"># el</span>

<span class="hljs-comment"># 常用方法</span>
<span class="hljs-built_in">print</span>(s.upper())          <span class="hljs-comment"># HELLO</span>
<span class="hljs-built_in">print</span>(s.replace(<span class="hljs-string">"l"</span>, <span class="hljs-string">"L"</span>))  <span class="hljs-comment"># heLLo</span>
<span class="hljs-built_in">print</span>(<span class="hljs-built_in">len</span>(s))             <span class="hljs-comment"># 5</span>
<span class="hljs-built_in">print</span>(<span class="hljs-string">", "</span>.join([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>]))  <span class="hljs-comment"># a, b, c</span>
</code></pre>
<h3>4.3 列表（list）</h3>
<p>列表是有序、可变的容器：</p>
<pre><code class="hljs language-python">fruits = [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>, <span class="hljs-string">"cherry"</span>]

fruits.append(<span class="hljs-string">"orange"</span>)      <span class="hljs-comment"># 末尾追加</span>
fruits.insert(<span class="hljs-number">0</span>, <span class="hljs-string">"grape"</span>)    <span class="hljs-comment"># 指定位置插入</span>
fruits.remove(<span class="hljs-string">"apple"</span>)       <span class="hljs-comment"># 按值删除</span>
popped = fruits.pop()        <span class="hljs-comment"># 弹出末尾元素</span>

<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">0</span>])             <span class="hljs-comment"># 索引</span>
<span class="hljs-built_in">print</span>(fruits[-<span class="hljs-number">1</span>])            <span class="hljs-comment"># 最后一个</span>
<span class="hljs-built_in">print</span>(fruits[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>])           <span class="hljs-comment"># 切片</span>

<span class="hljs-comment"># 列表推导式（非常常用）</span>
squares = [x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>)]
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]</span>
</code></pre>
<h3>4.4 元组（tuple）</h3>
<p>元组是不可变列表，常用于返回多个值：</p>
<pre><code class="hljs language-python">point = (<span class="hljs-number">3</span>, <span class="hljs-number">4</span>)
x, y = point          <span class="hljs-comment"># 解包</span>
<span class="hljs-built_in">print</span>(x, y)           <span class="hljs-comment"># 3 4</span>

<span class="hljs-keyword">def</span> <span class="hljs-title function_">minmax</span>(<span class="hljs-params">nums</span>):
    <span class="hljs-keyword">return</span> <span class="hljs-built_in">min</span>(nums), <span class="hljs-built_in">max</span>(nums)

lo, hi = minmax([<span class="hljs-number">3</span>, <span class="hljs-number">1</span>, <span class="hljs-number">4</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>])
<span class="hljs-built_in">print</span>(lo, hi)         <span class="hljs-comment"># 1 5</span>
</code></pre>
<h3>4.5 字典（dict）</h3>
<p>字典是键值对容器，查询效率高：</p>
<pre><code class="hljs language-python">person = {<span class="hljs-string">"name"</span>: <span class="hljs-string">"zorrooz"</span>, <span class="hljs-string">"age"</span>: <span class="hljs-number">25</span>}

<span class="hljs-built_in">print</span>(person[<span class="hljs-string">"name"</span>])            <span class="hljs-comment"># 取值</span>
<span class="hljs-built_in">print</span>(person.get(<span class="hljs-string">"city"</span>, <span class="hljs-string">"未知"</span>))  <span class="hljs-comment"># 安全的取值，带默认值</span>
person[<span class="hljs-string">"city"</span>] = <span class="hljs-string">"Lanzhou"</span>       <span class="hljs-comment"># 新增/修改</span>

<span class="hljs-keyword">for</span> key, value <span class="hljs-keyword">in</span> person.items():
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{key}</span>: <span class="hljs-subst">{value}</span>"</span>)

<span class="hljs-comment"># 字典推导式</span>
squares = {x: x ** <span class="hljs-number">2</span> <span class="hljs-keyword">for</span> x <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>)}
<span class="hljs-built_in">print</span>(squares)  <span class="hljs-comment"># {0: 0, 1: 1, 2: 4, 3: 9, 4: 16}</span>
</code></pre>
<h3>4.6 集合（set）</h3>
<p>集合是无序、去重的容器：</p>
<pre><code class="hljs language-python">a = {<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>, <span class="hljs-number">3</span>}
<span class="hljs-built_in">print</span>(a)              <span class="hljs-comment"># {1, 2, 3}，自动去重</span>

b = {<span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>}
<span class="hljs-built_in">print</span>(a &#x26; b)          <span class="hljs-comment"># 交集 {2, 3}</span>
<span class="hljs-built_in">print</span>(a | b)          <span class="hljs-comment"># 并集 {1, 2, 3, 4}</span>
<span class="hljs-built_in">print</span>(a - b)          <span class="hljs-comment"># 差集 {1}</span>
</code></pre>
<h2>5. 流程控制</h2>
<h3>5.1 条件判断</h3>
<pre><code class="hljs language-python">score = <span class="hljs-number">85</span>

<span class="hljs-keyword">if</span> score >= <span class="hljs-number">90</span>:
    grade = <span class="hljs-string">"A"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">80</span>:
    grade = <span class="hljs-string">"B"</span>
<span class="hljs-keyword">elif</span> score >= <span class="hljs-number">70</span>:
    grade = <span class="hljs-string">"C"</span>
<span class="hljs-keyword">else</span>:
    grade = <span class="hljs-string">"D"</span>

<span class="hljs-built_in">print</span>(<span class="hljs-string">f"等级：<span class="hljs-subst">{grade}</span>"</span>)

<span class="hljs-comment"># 三元表达式</span>
status = <span class="hljs-string">"pass"</span> <span class="hljs-keyword">if</span> score >= <span class="hljs-number">60</span> <span class="hljs-keyword">else</span> <span class="hljs-string">"fail"</span>
</code></pre>
<h3>5.2 for 循环</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 遍历可迭代对象</span>
<span class="hljs-keyword">for</span> fruit <span class="hljs-keyword">in</span> [<span class="hljs-string">"apple"</span>, <span class="hljs-string">"banana"</span>]:
    <span class="hljs-built_in">print</span>(fruit)

<span class="hljs-comment"># range 生成数列</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):        <span class="hljs-comment"># 0 到 4</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">1</span>, <span class="hljs-number">10</span>, <span class="hljs-number">2</span>):  <span class="hljs-comment"># 1, 3, 5, 7, 9</span>
    <span class="hljs-built_in">print</span>(i)

<span class="hljs-comment"># enumerate 同时取索引和值</span>
<span class="hljs-keyword">for</span> idx, fruit <span class="hljs-keyword">in</span> <span class="hljs-built_in">enumerate</span>([<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>]):
    <span class="hljs-built_in">print</span>(idx, fruit)
</code></pre>
<h3>5.3 while 循环</h3>
<pre><code class="hljs language-python">n = <span class="hljs-number">0</span>
<span class="hljs-keyword">while</span> n &#x3C; <span class="hljs-number">5</span>:
    <span class="hljs-built_in">print</span>(n)
    n += <span class="hljs-number">1</span>          <span class="hljs-comment"># 注意：忘记递增会死循环</span>
</code></pre>
<h3>5.4 break / continue / else</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># break：提前终止循环</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">10</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">3</span>:
        <span class="hljs-keyword">break</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># 输出 0 1 2</span>

<span class="hljs-comment"># continue：跳过本次迭代</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">5</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">2</span>:
        <span class="hljs-keyword">continue</span>
    <span class="hljs-built_in">print</span>(i)        <span class="hljs-comment"># 输出 0 1 3 4</span>

<span class="hljs-comment"># for-else：循环正常结束（未被 break）时执行</span>
<span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">3</span>):
    <span class="hljs-keyword">if</span> i == <span class="hljs-number">99</span>:
        <span class="hljs-keyword">break</span>
<span class="hljs-keyword">else</span>:
    <span class="hljs-built_in">print</span>(<span class="hljs-string">"循环正常结束"</span>)
</code></pre>
<h2>6. 综合练习</h2>
<p>统计一段文本中每个单词的出现次数：</p>
<pre><code class="hljs language-python">text = <span class="hljs-string">"python is fun and python is powerful"</span>

words = text.split()
counter = {}

<span class="hljs-keyword">for</span> word <span class="hljs-keyword">in</span> words:
    counter[word] = counter.get(word, <span class="hljs-number">0</span>) + <span class="hljs-number">1</span>

<span class="hljs-keyword">for</span> word, count <span class="hljs-keyword">in</span> <span class="hljs-built_in">sorted</span>(counter.items(), key=<span class="hljs-keyword">lambda</span> x: -x[<span class="hljs-number">1</span>]):
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{word}</span>: <span class="hljs-subst">{count}</span>"</span>)
</code></pre>
<p>输出：</p>
<pre><code>python: 2
is: 2
fun: 1
and: 1
powerful: 1
</code></pre>
<h2>7. 小结</h2>
<ul>
<li>用 <code>venv</code> 隔离项目依赖，用 <code>pip freeze</code> 记录依赖</li>
<li>内置类型：<code>int</code> / <code>float</code> / <code>str</code> / <code>list</code> / <code>tuple</code> / <code>dict</code> / <code>set</code></li>
<li>列表与字典推导式、f-string 是日常高频语法</li>
<li>流程控制：<code>if/elif/else</code>、<code>for</code>、<code>while</code>、<code>break</code> / <code>continue</code></li>
</ul>
<p>下一篇将介绍函数、类与模块，进入真正的工程化编程。</p>`,P_=`<h1>Python Data Processing in Practice: File IO, Regex, and NumPy/Pandas</h1>
<p>Data processing is the most core application scenario of Python. This article explains file I/O, regular expressions, and the two pillars of scientific computing, NumPy and Pandas, with examples tailored to common tasks in bioinformatics.</p>
<h2>1. File I/O</h2>
<h3>1.1 Text Files</h3>
<p>It is recommended to use the <code>with</code> statement to automatically manage file closure:</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># Write</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"w"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    f.write(<span class="hljs-string">"First line\\n"</span>)
    f.write(<span class="hljs-string">"Second line\\n"</span>)

<span class="hljs-comment"># Read all</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    content = f.read()
<span class="hljs-built_in">print</span>(content)

<span class="hljs-comment"># Read line by line (recommended for large files)</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
        <span class="hljs-built_in">print</span>(line.strip())   <span class="hljs-comment"># strip removes the newline</span>
</code></pre>
<p>Mode description: <code>r</code> read, <code>w</code> write (overwrite), <code>a</code> append, <code>rb</code>/<code>wb</code> binary.</p>
<h3>1.2 Parsing FASTA Files</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">read_fasta</span>(<span class="hljs-params">path</span>):
    <span class="hljs-string">"""Read a FASTA file and return a {id: sequence} dictionary."""</span>
    sequences = {}
    current_id = <span class="hljs-literal">None</span>
    <span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(path, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
        <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
            line = line.strip()
            <span class="hljs-keyword">if</span> line.startswith(<span class="hljs-string">">"</span>):
                current_id = line[<span class="hljs-number">1</span>:].split()[<span class="hljs-number">0</span>]   <span class="hljs-comment"># take the first field as the id</span>
                sequences[current_id] = <span class="hljs-string">""</span>
            <span class="hljs-keyword">elif</span> current_id <span class="hljs-keyword">is</span> <span class="hljs-keyword">not</span> <span class="hljs-literal">None</span>:
                sequences[current_id] += line
    <span class="hljs-keyword">return</span> sequences

seqs = read_fasta(<span class="hljs-string">"example.fasta"</span>)
<span class="hljs-keyword">for</span> seq_id, seq <span class="hljs-keyword">in</span> seqs.items():
    <span class="hljs-built_in">print</span>(seq_id, <span class="hljs-built_in">len</span>(seq))
</code></pre>
<h3>1.3 CSV Files</h3>
<p>Standard library <code>csv</code> module:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> csv

<span class="hljs-comment"># Read</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"data.csv"</span>, <span class="hljs-string">"r"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    reader = csv.DictReader(f)          <span class="hljs-comment"># read as dict by header</span>
    rows = <span class="hljs-built_in">list</span>(reader)
    <span class="hljs-built_in">print</span>(rows[<span class="hljs-number">0</span>][<span class="hljs-string">"gene"</span>], rows[<span class="hljs-number">0</span>][<span class="hljs-string">"expression"</span>])

<span class="hljs-comment"># Write</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"out.csv"</span>, <span class="hljs-string">"w"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    writer = csv.writer(f)
    writer.writerow([<span class="hljs-string">"gene"</span>, <span class="hljs-string">"value"</span>])
    writer.writerows([[<span class="hljs-string">"TP53"</span>, <span class="hljs-number">12.3</span>], [<span class="hljs-string">"BRCA1"</span>, <span class="hljs-number">8.9</span>]])
</code></pre>
<blockquote>
<p>For big data scenarios, Pandas' <code>read_csv</code> is recommended (see Section 4).</p>
</blockquote>
<h2>2. Regular Expressions</h2>
<p>Regular expressions are used for text pattern matching. Common functions in the <code>re</code> module:</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> re

text = <span class="hljs-string">"The protein P12345 has 500 amino acids, code: A1B2C3"</span>

<span class="hljs-comment"># match: match from the beginning</span>
m = re.<span class="hljs-keyword">match</span>(<span class="hljs-string">r"The"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># The</span>

<span class="hljs-comment"># search: search for the first match</span>
m = re.search(<span class="hljs-string">r"\\d+"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># 12345 (the first digit string)</span>

<span class="hljs-comment"># findall: return all matches</span>
<span class="hljs-built_in">print</span>(re.findall(<span class="hljs-string">r"[A-Z]\\d"</span>, text))   <span class="hljs-comment"># ['P1', 'A1', 'B2', 'C3']</span>

<span class="hljs-comment"># sub: substitute</span>
<span class="hljs-built_in">print</span>(re.sub(<span class="hljs-string">r"\\d+"</span>, <span class="hljs-string">"#"</span>, text))
</code></pre>
<h3>2.1 Common Metacharacters</h3>
<table>
<thead>
<tr>
<th>Pattern</th>
<th>Meaning</th>
<th>Example</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>.</code></td>
<td>Any character (except newline)</td>
<td><code>a.c</code> matches <code>abc</code></td>
</tr>
<tr>
<td><code>*</code></td>
<td>Previous character 0 or more times</td>
<td><code>ab*c</code> matches <code>ac</code>, <code>abc</code></td>
</tr>
<tr>
<td><code>+</code></td>
<td>Previous character 1 or more times</td>
<td><code>ab+c</code> matches <code>abc</code></td>
</tr>
<tr>
<td><code>?</code></td>
<td>Previous character 0 or 1 time</td>
<td><code>colou?r</code> matches color/colour</td>
</tr>
<tr>
<td><code>^</code> / <code>$</code></td>
<td>Start / end of line</td>
<td><code>^ATG</code> starts with ATG</td>
</tr>
<tr>
<td><code>[abc]</code></td>
<td>Character set</td>
<td><code>[ATGC]</code> any base</td>
</tr>
<tr>
<td><code>\\d</code> / <code>\\w</code> / <code>\\s</code></td>
<td>Digit / word character / whitespace</td>
<td></td>
</tr>
<tr>
<td><code>{n,m}</code></td>
<td>Repeat n to m times</td>
<td><code>\\d{3}</code> three digits</td>
</tr>
<tr>
<td><code>(…)</code></td>
<td>Group capture</td>
<td><code>(ATG){2}</code></td>
</tr>
</tbody>
</table>
<h3>2.2 Practical Information Extraction</h3>
<p>Extract coordinates and names from gene annotation text:</p>
<pre><code class="hljs language-python">line = <span class="hljs-string">'gene=TP53;gene_id=ENSG00000141510;chromosome=17;start=7668402;end=7687550'</span>

m = re.search(<span class="hljs-string">r"gene=(\\w+);.*?start=(\\d+);end=(\\d+)"</span>, line)
<span class="hljs-keyword">if</span> m:
    gene, start, end = m.groups()
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{gene}</span>: <span class="hljs-subst">{start}</span>-<span class="hljs-subst">{end}</span> (<span class="hljs-subst">{<span class="hljs-built_in">int</span>(end) - <span class="hljs-built_in">int</span>(start)}</span> bp)"</span>)
<span class="hljs-comment"># TP53: 7668402-7687550 (19148 bp)</span>
</code></pre>
<h3>2.3 Compiling Regular Expressions for Performance</h3>
<pre><code class="hljs language-python">motif = re.<span class="hljs-built_in">compile</span>(<span class="hljs-string">r"^ATG.*(TAA|TAG|TGA)$"</span>)   <span class="hljs-comment"># rough match for complete ORF</span>
<span class="hljs-keyword">for</span> seq <span class="hljs-keyword">in</span> [<span class="hljs-string">"ATGCCCTAA"</span>, <span class="hljs-string">"ATGGGG"</span>]:
    <span class="hljs-built_in">print</span>(seq, <span class="hljs-built_in">bool</span>(motif.<span class="hljs-keyword">match</span>(seq)))
</code></pre>
<h2>3. NumPy: Basics of Numerical Computing</h2>
<h3>3.1 Array Creation</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np

a = np.array([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>])
b = np.zeros((<span class="hljs-number">2</span>, <span class="hljs-number">3</span>))          <span class="hljs-comment"># 2x3 all-zero matrix</span>
c = np.ones(<span class="hljs-number">5</span>)
d = np.arange(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">0.2</span>)      <span class="hljs-comment"># [0.  0.2 0.4 0.6 0.8]</span>
e = np.linspace(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>)      <span class="hljs-comment"># split 0 to 1 into 5 equal parts</span>
f = np.random.rand(<span class="hljs-number">3</span>, <span class="hljs-number">3</span>)      <span class="hljs-comment"># uniformly distributed random numbers</span>
</code></pre>
<h3>3.2 Vectorized Operations (Core Advantage)</h3>
<pre><code class="hljs language-python">x = np.array([<span class="hljs-number">1.0</span>, <span class="hljs-number">2.0</span>, <span class="hljs-number">3.0</span>])
y = np.array([<span class="hljs-number">4.0</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">6.0</span>])

<span class="hljs-built_in">print</span>(x + y)          <span class="hljs-comment"># [5. 7. 9.] without loops</span>
<span class="hljs-built_in">print</span>(x * y)          <span class="hljs-comment"># [4. 10. 18.]</span>
<span class="hljs-built_in">print</span>(np.sqrt(x))     <span class="hljs-comment"># element-wise square root</span>
<span class="hljs-built_in">print</span>(x.mean(), x.std())   <span class="hljs-comment"># 2.0 0.816...</span>

<span class="hljs-comment"># Broadcasting: automatic expansion of different shapes</span>
matrix = np.array([[<span class="hljs-number">1</span>, <span class="hljs-number">2</span>], [<span class="hljs-number">3</span>, <span class="hljs-number">4</span>]])
<span class="hljs-built_in">print</span>(matrix * <span class="hljs-number">10</span>)            <span class="hljs-comment"># multiply each element by 10</span>
<span class="hljs-built_in">print</span>(matrix + np.array([<span class="hljs-number">100</span>, <span class="hljs-number">200</span>]))  <span class="hljs-comment"># row broadcast</span>
</code></pre>
<h3>3.3 Indexing and Slicing</h3>
<pre><code class="hljs language-python">arr = np.arange(<span class="hljs-number">10</span>)
<span class="hljs-built_in">print</span>(arr[<span class="hljs-number">2</span>:<span class="hljs-number">8</span>:<span class="hljs-number">2</span>])        <span class="hljs-comment"># [2 4 6], start:stop:step</span>
<span class="hljs-built_in">print</span>(arr[arr > <span class="hljs-number">5</span>])      <span class="hljs-comment"># boolean mask [6 7 8 9]</span>

m = np.random.rand(<span class="hljs-number">4</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">0</span>, :])           <span class="hljs-comment"># first row</span>
<span class="hljs-built_in">print</span>(m[:, <span class="hljs-number">1</span>])           <span class="hljs-comment"># second column</span>
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>, <span class="hljs-number">1</span>:<span class="hljs-number">3</span>])       <span class="hljs-comment"># submatrix</span>
</code></pre>
<h3>3.4 Bioinformatics Example: Normalizing an Expression Matrix</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Simulate a gene expression matrix: 4 genes × 5 samples</span>
expr = np.array([
    [<span class="hljs-number">12.3</span>, <span class="hljs-number">15.1</span>, <span class="hljs-number">9.8</span>, <span class="hljs-number">20.4</span>, <span class="hljs-number">11.2</span>],
    [<span class="hljs-number">3.2</span>, <span class="hljs-number">4.1</span>, <span class="hljs-number">2.9</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">3.6</span>],
    [<span class="hljs-number">80.5</span>, <span class="hljs-number">92.3</span>, <span class="hljs-number">75.1</span>, <span class="hljs-number">100.2</span>, <span class="hljs-number">88.7</span>],
    [<span class="hljs-number">45.6</span>, <span class="hljs-number">41.2</span>, <span class="hljs-number">48.9</span>, <span class="hljs-number">39.8</span>, <span class="hljs-number">50.1</span>],
])

<span class="hljs-comment"># Perform Z-score normalization by row (gene): (x - mean) / std</span>
means = expr.mean(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
stds = expr.std(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
z = (expr - means) / stds

<span class="hljs-built_in">print</span>(np.<span class="hljs-built_in">round</span>(z, <span class="hljs-number">2</span>))
</code></pre>
<h2>4. Pandas: Tabular Data Processing</h2>
<h3>4.1 Series and DataFrame</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> pandas <span class="hljs-keyword">as</span> pd

<span class="hljs-comment"># Series: one-dimensional labeled array</span>
s = pd.Series([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>], index=[<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>])

<span class="hljs-comment"># DataFrame: two-dimensional table</span>
df = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"BRCA1"</span>, <span class="hljs-string">"EGFR"</span>, <span class="hljs-string">"MYC"</span>],
    <span class="hljs-string">"expression"</span>: [<span class="hljs-number">12.3</span>, <span class="hljs-number">8.9</span>, <span class="hljs-number">45.1</span>, <span class="hljs-number">33.2</span>],
    <span class="hljs-string">"chromosome"</span>: [<span class="hljs-number">17</span>, <span class="hljs-number">17</span>, <span class="hljs-number">7</span>, <span class="hljs-number">8</span>],
})
<span class="hljs-built_in">print</span>(df)
</code></pre>
<h3>4.2 Reading and Writing</h3>
<pre><code class="hljs language-python">df = pd.read_csv(<span class="hljs-string">"expression.csv"</span>)
df = pd.read_csv(<span class="hljs-string">"data.tsv"</span>, sep=<span class="hljs-string">"\\t"</span>)
df = pd.read_excel(<span class="hljs-string">"data.xlsx"</span>, sheet_name=<span class="hljs-string">"Sheet1"</span>)
df = pd.read_json(<span class="hljs-string">"data.json"</span>)

df.to_csv(<span class="hljs-string">"out.csv"</span>, index=<span class="hljs-literal">False</span>)
</code></pre>
<h3>4.3 Viewing and Filtering</h3>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(df.head())          <span class="hljs-comment"># first 5 rows</span>
<span class="hljs-built_in">print</span>(df.info())          <span class="hljs-comment"># column types and missing values</span>
<span class="hljs-built_in">print</span>(df.describe())      <span class="hljs-comment"># statistical summary of numeric columns</span>

<span class="hljs-comment"># Column selection</span>
<span class="hljs-built_in">print</span>(df[<span class="hljs-string">"gene"</span>])
<span class="hljs-built_in">print</span>(df[[<span class="hljs-string">"gene"</span>, <span class="hljs-string">"expression"</span>]])

<span class="hljs-comment"># Row filtering (boolean indexing)</span>
high = df[df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>]
<span class="hljs-built_in">print</span>(high)

<span class="hljs-comment"># Multiple conditions</span>
sel = df[(df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>) &#x26; (df[<span class="hljs-string">"chromosome"</span>] == <span class="hljs-number">17</span>)]
<span class="hljs-built_in">print</span>(sel)

<span class="hljs-comment"># isin filtering</span>
sel = df[df[<span class="hljs-string">"gene"</span>].isin([<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"MYC"</span>])]
</code></pre>
<h3>4.4 Grouping and Aggregation</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Group by chromosome and compute expression mean and count</span>
grouped = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].agg([<span class="hljs-string">"mean"</span>, <span class="hljs-string">"count"</span>, <span class="hljs-string">"max"</span>])
<span class="hljs-built_in">print</span>(grouped)

<span class="hljs-comment"># Transform by group (e.g., within-group standardization)</span>
df[<span class="hljs-string">"zscore"</span>] = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].transform(
    <span class="hljs-keyword">lambda</span> x: (x - x.mean()) / x.std()
)
</code></pre>
<h3>4.5 Handling Missing Values</h3>
<pre><code class="hljs language-python">df = pd.DataFrame({<span class="hljs-string">"a"</span>: [<span class="hljs-number">1</span>, <span class="hljs-literal">None</span>, <span class="hljs-number">3</span>], <span class="hljs-string">"b"</span>: [<span class="hljs-number">4</span>, <span class="hljs-number">5</span>, <span class="hljs-literal">None</span>]})

<span class="hljs-built_in">print</span>(df.isna().<span class="hljs-built_in">sum</span>())        <span class="hljs-comment"># missing value count per column</span>
<span class="hljs-built_in">print</span>(df.dropna())            <span class="hljs-comment"># drop rows with missing values</span>
<span class="hljs-built_in">print</span>(df.fillna(<span class="hljs-number">0</span>))           <span class="hljs-comment"># fill with 0</span>
<span class="hljs-built_in">print</span>(df.fillna(df.mean()))   <span class="hljs-comment"># fill with mean</span>
</code></pre>
<h3>4.6 Practical Example: Filtering Differentially Expressed Genes</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Simulate: gene + treated/control expression + p-value</span>
results = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">f"GENE<span class="hljs-subst">{i}</span>"</span> <span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">100</span>)],
    <span class="hljs-string">"log2fc"</span>: np.random.normal(<span class="hljs-number">0</span>, <span class="hljs-number">1.5</span>, <span class="hljs-number">100</span>),
    <span class="hljs-string">"pvalue"</span>: np.random.rand(<span class="hljs-number">100</span>),
})

<span class="hljs-comment"># Filter significant genes: |log2FC| > 1 and p &#x3C; 0.05</span>
sig = results[(results[<span class="hljs-string">"log2fc"</span>].<span class="hljs-built_in">abs</span>() > <span class="hljs-number">1</span>) &#x26; (results[<span class="hljs-string">"pvalue"</span>] &#x3C; <span class="hljs-number">0.05</span>)]

<span class="hljs-comment"># Sort by log2FC</span>
sig = sig.sort_values(<span class="hljs-string">"log2fc"</span>, ascending=<span class="hljs-literal">False</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"Significant gene count: <span class="hljs-subst">{<span class="hljs-built_in">len</span>(sig)}</span>"</span>)
<span class="hljs-built_in">print</span>(sig.head())
</code></pre>
<h2>5. Summary</h2>
<ul>
<li>Files: <code>with open</code> + line-by-line reading for large files; use <code>csv</code> module or Pandas for CSV</li>
<li>Regex: <code>re.search</code> / <code>findall</code> / <code>sub</code>, start with small patterns and combine</li>
<li>NumPy: vectorized operations + broadcasting, avoid Python loops</li>
<li>Pandas: <code>read_csv</code> → filtering/grouping/missing value handling → export, the standard workflow for tabular data</li>
</ul>
<p>At this point, the Python tutorial trilogy is complete: Beginner → Intermediate → Data Processing.</p>`,S_=`<h1>Python 数据处理实战：文件 IO、正则与 NumPy/Pandas</h1>
<p>数据处理是 Python 最核心的应用场景。本文讲解文件读写、正则表达式，以及科学计算的两大支柱 NumPy 与 Pandas，示例贴合生物信息学常见任务。</p>
<h2>1. 文件读写</h2>
<h3>1.1 文本文件</h3>
<p>推荐使用 <code>with</code> 语句，自动管理文件关闭：</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># 写入</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"w"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    f.write(<span class="hljs-string">"第一行\\n"</span>)
    f.write(<span class="hljs-string">"第二行\\n"</span>)

<span class="hljs-comment"># 读取全部</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    content = f.read()
<span class="hljs-built_in">print</span>(content)

<span class="hljs-comment"># 逐行读取（大文件推荐）</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"output.txt"</span>, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
        <span class="hljs-built_in">print</span>(line.strip())   <span class="hljs-comment"># strip 去掉换行符</span>
</code></pre>
<p>模式说明：<code>r</code> 读、<code>w</code> 写（覆盖）、<code>a</code> 追加、<code>rb</code>/<code>wb</code> 二进制。</p>
<h3>1.2 解析 FASTA 文件</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">def</span> <span class="hljs-title function_">read_fasta</span>(<span class="hljs-params">path</span>):
    <span class="hljs-string">"""读取 FASTA 文件，返回 {id: 序列} 字典。"""</span>
    sequences = {}
    current_id = <span class="hljs-literal">None</span>
    <span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(path, <span class="hljs-string">"r"</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
        <span class="hljs-keyword">for</span> line <span class="hljs-keyword">in</span> f:
            line = line.strip()
            <span class="hljs-keyword">if</span> line.startswith(<span class="hljs-string">">"</span>):
                current_id = line[<span class="hljs-number">1</span>:].split()[<span class="hljs-number">0</span>]   <span class="hljs-comment"># 取第一个字段作为 id</span>
                sequences[current_id] = <span class="hljs-string">""</span>
            <span class="hljs-keyword">elif</span> current_id <span class="hljs-keyword">is</span> <span class="hljs-keyword">not</span> <span class="hljs-literal">None</span>:
                sequences[current_id] += line
    <span class="hljs-keyword">return</span> sequences

seqs = read_fasta(<span class="hljs-string">"example.fasta"</span>)
<span class="hljs-keyword">for</span> seq_id, seq <span class="hljs-keyword">in</span> seqs.items():
    <span class="hljs-built_in">print</span>(seq_id, <span class="hljs-built_in">len</span>(seq))
</code></pre>
<h3>1.3 CSV 文件</h3>
<p>标准库 <code>csv</code> 模块：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> csv

<span class="hljs-comment"># 读取</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"data.csv"</span>, <span class="hljs-string">"r"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    reader = csv.DictReader(f)          <span class="hljs-comment"># 按表头读取为字典</span>
    rows = <span class="hljs-built_in">list</span>(reader)
    <span class="hljs-built_in">print</span>(rows[<span class="hljs-number">0</span>][<span class="hljs-string">"gene"</span>], rows[<span class="hljs-number">0</span>][<span class="hljs-string">"expression"</span>])

<span class="hljs-comment"># 写入</span>
<span class="hljs-keyword">with</span> <span class="hljs-built_in">open</span>(<span class="hljs-string">"out.csv"</span>, <span class="hljs-string">"w"</span>, newline=<span class="hljs-string">""</span>, encoding=<span class="hljs-string">"utf-8"</span>) <span class="hljs-keyword">as</span> f:
    writer = csv.writer(f)
    writer.writerow([<span class="hljs-string">"gene"</span>, <span class="hljs-string">"value"</span>])
    writer.writerows([[<span class="hljs-string">"TP53"</span>, <span class="hljs-number">12.3</span>], [<span class="hljs-string">"BRCA1"</span>, <span class="hljs-number">8.9</span>]])
</code></pre>
<blockquote>
<p>大数据场景推荐 Pandas 的 <code>read_csv</code>（见第 4 节）。</p>
</blockquote>
<h2>2. 正则表达式</h2>
<p>正则用于文本模式匹配。<code>re</code> 模块常用函数：</p>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> re

text = <span class="hljs-string">"The protein P12345 has 500 amino acids, code: A1B2C3"</span>

<span class="hljs-comment"># match：从开头匹配</span>
m = re.<span class="hljs-keyword">match</span>(<span class="hljs-string">r"The"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># The</span>

<span class="hljs-comment"># search：搜索第一个匹配</span>
m = re.search(<span class="hljs-string">r"\\d+"</span>, text)
<span class="hljs-built_in">print</span>(m.group())          <span class="hljs-comment"># 12345（第一个数字串）</span>

<span class="hljs-comment"># findall：返回所有匹配</span>
<span class="hljs-built_in">print</span>(re.findall(<span class="hljs-string">r"[A-Z]\\d"</span>, text))   <span class="hljs-comment"># ['P1', 'A1', 'B2', 'C3']</span>

<span class="hljs-comment"># sub：替换</span>
<span class="hljs-built_in">print</span>(re.sub(<span class="hljs-string">r"\\d+"</span>, <span class="hljs-string">"#"</span>, text))
</code></pre>
<h3>2.1 常用元字符</h3>
<table>
<thead>
<tr>
<th>模式</th>
<th>含义</th>
<th>示例</th>
</tr>
</thead>
<tbody>
<tr>
<td><code>.</code></td>
<td>任意字符（除换行）</td>
<td><code>a.c</code> 匹配 <code>abc</code></td>
</tr>
<tr>
<td><code>*</code></td>
<td>前一个字符 0 次或多次</td>
<td><code>ab*c</code> 匹配 <code>ac</code>、<code>abc</code></td>
</tr>
<tr>
<td><code>+</code></td>
<td>前一个字符 1 次或多次</td>
<td><code>ab+c</code> 匹配 <code>abc</code></td>
</tr>
<tr>
<td><code>?</code></td>
<td>前一个字符 0 次或 1 次</td>
<td><code>colou?r</code> 匹配 color/colour</td>
</tr>
<tr>
<td><code>^</code> / <code>$</code></td>
<td>行首 / 行尾</td>
<td><code>^ATG</code> 以 ATG 开头</td>
</tr>
<tr>
<td><code>[abc]</code></td>
<td>字符集</td>
<td><code>[ATGC]</code> 任意碱基</td>
</tr>
<tr>
<td><code>\\d</code> / <code>\\w</code> / <code>\\s</code></td>
<td>数字 / 单词字符 / 空白</td>
<td></td>
</tr>
<tr>
<td><code>{n,m}</code></td>
<td>重复 n 到 m 次</td>
<td><code>\\d{3}</code> 三个数字</td>
</tr>
<tr>
<td><code>(…)</code></td>
<td>分组捕获</td>
<td><code>(ATG){2}</code></td>
</tr>
</tbody>
</table>
<h3>2.2 提取信息实战</h3>
<p>从基因注释文本提取坐标与名称：</p>
<pre><code class="hljs language-python">line = <span class="hljs-string">'gene=TP53;gene_id=ENSG00000141510;chromosome=17;start=7668402;end=7687550'</span>

m = re.search(<span class="hljs-string">r"gene=(\\w+);.*?start=(\\d+);end=(\\d+)"</span>, line)
<span class="hljs-keyword">if</span> m:
    gene, start, end = m.groups()
    <span class="hljs-built_in">print</span>(<span class="hljs-string">f"<span class="hljs-subst">{gene}</span>: <span class="hljs-subst">{start}</span>-<span class="hljs-subst">{end}</span> (<span class="hljs-subst">{<span class="hljs-built_in">int</span>(end) - <span class="hljs-built_in">int</span>(start)}</span> bp)"</span>)
<span class="hljs-comment"># TP53: 7668402-7687550 (19148 bp)</span>
</code></pre>
<h3>2.3 编译正则提升性能</h3>
<pre><code class="hljs language-python">motif = re.<span class="hljs-built_in">compile</span>(<span class="hljs-string">r"^ATG.*(TAA|TAG|TGA)$"</span>)   <span class="hljs-comment"># 完整 ORF 粗略匹配</span>
<span class="hljs-keyword">for</span> seq <span class="hljs-keyword">in</span> [<span class="hljs-string">"ATGCCCTAA"</span>, <span class="hljs-string">"ATGGGG"</span>]:
    <span class="hljs-built_in">print</span>(seq, <span class="hljs-built_in">bool</span>(motif.<span class="hljs-keyword">match</span>(seq)))
</code></pre>
<h2>3. NumPy：数值计算基础</h2>
<h3>3.1 数组创建</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> numpy <span class="hljs-keyword">as</span> np

a = np.array([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>, <span class="hljs-number">4</span>])
b = np.zeros((<span class="hljs-number">2</span>, <span class="hljs-number">3</span>))          <span class="hljs-comment"># 2x3 全零矩阵</span>
c = np.ones(<span class="hljs-number">5</span>)
d = np.arange(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">0.2</span>)      <span class="hljs-comment"># [0.  0.2 0.4 0.6 0.8]</span>
e = np.linspace(<span class="hljs-number">0</span>, <span class="hljs-number">1</span>, <span class="hljs-number">5</span>)      <span class="hljs-comment"># 0 到 1 等分 5 份</span>
f = np.random.rand(<span class="hljs-number">3</span>, <span class="hljs-number">3</span>)      <span class="hljs-comment"># 均匀分布随机数</span>
</code></pre>
<h3>3.2 向量化运算（核心优势）</h3>
<pre><code class="hljs language-python">x = np.array([<span class="hljs-number">1.0</span>, <span class="hljs-number">2.0</span>, <span class="hljs-number">3.0</span>])
y = np.array([<span class="hljs-number">4.0</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">6.0</span>])

<span class="hljs-built_in">print</span>(x + y)          <span class="hljs-comment"># [5. 7. 9.]，无需循环</span>
<span class="hljs-built_in">print</span>(x * y)          <span class="hljs-comment"># [4. 10. 18.]</span>
<span class="hljs-built_in">print</span>(np.sqrt(x))     <span class="hljs-comment"># 逐元素开方</span>
<span class="hljs-built_in">print</span>(x.mean(), x.std())   <span class="hljs-comment"># 2.0 0.816...</span>

<span class="hljs-comment"># 广播：不同形状自动扩展</span>
matrix = np.array([[<span class="hljs-number">1</span>, <span class="hljs-number">2</span>], [<span class="hljs-number">3</span>, <span class="hljs-number">4</span>]])
<span class="hljs-built_in">print</span>(matrix * <span class="hljs-number">10</span>)            <span class="hljs-comment"># 每个元素乘 10</span>
<span class="hljs-built_in">print</span>(matrix + np.array([<span class="hljs-number">100</span>, <span class="hljs-number">200</span>]))  <span class="hljs-comment"># 行广播</span>
</code></pre>
<h3>3.3 索引与切片</h3>
<pre><code class="hljs language-python">arr = np.arange(<span class="hljs-number">10</span>)
<span class="hljs-built_in">print</span>(arr[<span class="hljs-number">2</span>:<span class="hljs-number">8</span>:<span class="hljs-number">2</span>])        <span class="hljs-comment"># [2 4 6]，起始:结束:步长</span>
<span class="hljs-built_in">print</span>(arr[arr > <span class="hljs-number">5</span>])      <span class="hljs-comment"># 布尔掩码 [6 7 8 9]</span>

m = np.random.rand(<span class="hljs-number">4</span>, <span class="hljs-number">4</span>)
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">0</span>, :])           <span class="hljs-comment"># 第一行</span>
<span class="hljs-built_in">print</span>(m[:, <span class="hljs-number">1</span>])           <span class="hljs-comment"># 第二列</span>
<span class="hljs-built_in">print</span>(m[<span class="hljs-number">1</span>:<span class="hljs-number">3</span>, <span class="hljs-number">1</span>:<span class="hljs-number">3</span>])       <span class="hljs-comment"># 子矩阵</span>
</code></pre>
<h3>3.4 生物信息学示例：归一化表达矩阵</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 模拟基因表达矩阵：4 基因 × 5 样本</span>
expr = np.array([
    [<span class="hljs-number">12.3</span>, <span class="hljs-number">15.1</span>, <span class="hljs-number">9.8</span>, <span class="hljs-number">20.4</span>, <span class="hljs-number">11.2</span>],
    [<span class="hljs-number">3.2</span>, <span class="hljs-number">4.1</span>, <span class="hljs-number">2.9</span>, <span class="hljs-number">5.0</span>, <span class="hljs-number">3.6</span>],
    [<span class="hljs-number">80.5</span>, <span class="hljs-number">92.3</span>, <span class="hljs-number">75.1</span>, <span class="hljs-number">100.2</span>, <span class="hljs-number">88.7</span>],
    [<span class="hljs-number">45.6</span>, <span class="hljs-number">41.2</span>, <span class="hljs-number">48.9</span>, <span class="hljs-number">39.8</span>, <span class="hljs-number">50.1</span>],
])

<span class="hljs-comment"># 按行（基因）做 Z-score 标准化：(x - mean) / std</span>
means = expr.mean(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
stds = expr.std(axis=<span class="hljs-number">1</span>, keepdims=<span class="hljs-literal">True</span>)
z = (expr - means) / stds

<span class="hljs-built_in">print</span>(np.<span class="hljs-built_in">round</span>(z, <span class="hljs-number">2</span>))
</code></pre>
<h2>4. Pandas：表格数据处理</h2>
<h3>4.1 Series 与 DataFrame</h3>
<pre><code class="hljs language-python"><span class="hljs-keyword">import</span> pandas <span class="hljs-keyword">as</span> pd

<span class="hljs-comment"># Series：一维带标签数组</span>
s = pd.Series([<span class="hljs-number">1</span>, <span class="hljs-number">2</span>, <span class="hljs-number">3</span>], index=[<span class="hljs-string">"a"</span>, <span class="hljs-string">"b"</span>, <span class="hljs-string">"c"</span>])

<span class="hljs-comment"># DataFrame：二维表格</span>
df = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"BRCA1"</span>, <span class="hljs-string">"EGFR"</span>, <span class="hljs-string">"MYC"</span>],
    <span class="hljs-string">"expression"</span>: [<span class="hljs-number">12.3</span>, <span class="hljs-number">8.9</span>, <span class="hljs-number">45.1</span>, <span class="hljs-number">33.2</span>],
    <span class="hljs-string">"chromosome"</span>: [<span class="hljs-number">17</span>, <span class="hljs-number">17</span>, <span class="hljs-number">7</span>, <span class="hljs-number">8</span>],
})
<span class="hljs-built_in">print</span>(df)
</code></pre>
<h3>4.2 读取与写入</h3>
<pre><code class="hljs language-python">df = pd.read_csv(<span class="hljs-string">"expression.csv"</span>)
df = pd.read_csv(<span class="hljs-string">"data.tsv"</span>, sep=<span class="hljs-string">"\\t"</span>)
df = pd.read_excel(<span class="hljs-string">"data.xlsx"</span>, sheet_name=<span class="hljs-string">"Sheet1"</span>)
df = pd.read_json(<span class="hljs-string">"data.json"</span>)

df.to_csv(<span class="hljs-string">"out.csv"</span>, index=<span class="hljs-literal">False</span>)
</code></pre>
<h3>4.3 查看与筛选</h3>
<pre><code class="hljs language-python"><span class="hljs-built_in">print</span>(df.head())          <span class="hljs-comment"># 前 5 行</span>
<span class="hljs-built_in">print</span>(df.info())          <span class="hljs-comment"># 列类型与缺失值</span>
<span class="hljs-built_in">print</span>(df.describe())      <span class="hljs-comment"># 数值列统计摘要</span>

<span class="hljs-comment"># 列选择</span>
<span class="hljs-built_in">print</span>(df[<span class="hljs-string">"gene"</span>])
<span class="hljs-built_in">print</span>(df[[<span class="hljs-string">"gene"</span>, <span class="hljs-string">"expression"</span>]])

<span class="hljs-comment"># 行筛选（布尔索引）</span>
high = df[df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>]
<span class="hljs-built_in">print</span>(high)

<span class="hljs-comment"># 多条件</span>
sel = df[(df[<span class="hljs-string">"expression"</span>] > <span class="hljs-number">10</span>) &#x26; (df[<span class="hljs-string">"chromosome"</span>] == <span class="hljs-number">17</span>)]
<span class="hljs-built_in">print</span>(sel)

<span class="hljs-comment"># isin 筛选</span>
sel = df[df[<span class="hljs-string">"gene"</span>].isin([<span class="hljs-string">"TP53"</span>, <span class="hljs-string">"MYC"</span>])]
</code></pre>
<h3>4.4 分组与聚合</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 按染色体分组，统计表达均值与数量</span>
grouped = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].agg([<span class="hljs-string">"mean"</span>, <span class="hljs-string">"count"</span>, <span class="hljs-string">"max"</span>])
<span class="hljs-built_in">print</span>(grouped)

<span class="hljs-comment"># 按分组转换（如组内标准化）</span>
df[<span class="hljs-string">"zscore"</span>] = df.groupby(<span class="hljs-string">"chromosome"</span>)[<span class="hljs-string">"expression"</span>].transform(
    <span class="hljs-keyword">lambda</span> x: (x - x.mean()) / x.std()
)
</code></pre>
<h3>4.5 缺失值处理</h3>
<pre><code class="hljs language-python">df = pd.DataFrame({<span class="hljs-string">"a"</span>: [<span class="hljs-number">1</span>, <span class="hljs-literal">None</span>, <span class="hljs-number">3</span>], <span class="hljs-string">"b"</span>: [<span class="hljs-number">4</span>, <span class="hljs-number">5</span>, <span class="hljs-literal">None</span>]})

<span class="hljs-built_in">print</span>(df.isna().<span class="hljs-built_in">sum</span>())        <span class="hljs-comment"># 每列缺失数量</span>
<span class="hljs-built_in">print</span>(df.dropna())            <span class="hljs-comment"># 删除含缺失的行</span>
<span class="hljs-built_in">print</span>(df.fillna(<span class="hljs-number">0</span>))           <span class="hljs-comment"># 填充为 0</span>
<span class="hljs-built_in">print</span>(df.fillna(df.mean()))   <span class="hljs-comment"># 用均值填充</span>
</code></pre>
<h3>4.6 实战：差异表达基因筛选</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 模拟：基因 + 处理组/对照组表达 + p 值</span>
results = pd.DataFrame({
    <span class="hljs-string">"gene"</span>: [<span class="hljs-string">f"GENE<span class="hljs-subst">{i}</span>"</span> <span class="hljs-keyword">for</span> i <span class="hljs-keyword">in</span> <span class="hljs-built_in">range</span>(<span class="hljs-number">100</span>)],
    <span class="hljs-string">"log2fc"</span>: np.random.normal(<span class="hljs-number">0</span>, <span class="hljs-number">1.5</span>, <span class="hljs-number">100</span>),
    <span class="hljs-string">"pvalue"</span>: np.random.rand(<span class="hljs-number">100</span>),
})

<span class="hljs-comment"># 筛选显著差异基因：|log2FC| > 1 且 p &#x3C; 0.05</span>
sig = results[(results[<span class="hljs-string">"log2fc"</span>].<span class="hljs-built_in">abs</span>() > <span class="hljs-number">1</span>) &#x26; (results[<span class="hljs-string">"pvalue"</span>] &#x3C; <span class="hljs-number">0.05</span>)]

<span class="hljs-comment"># 按 log2FC 排序</span>
sig = sig.sort_values(<span class="hljs-string">"log2fc"</span>, ascending=<span class="hljs-literal">False</span>)
<span class="hljs-built_in">print</span>(<span class="hljs-string">f"显著基因数: <span class="hljs-subst">{<span class="hljs-built_in">len</span>(sig)}</span>"</span>)
<span class="hljs-built_in">print</span>(sig.head())
</code></pre>
<h2>5. 小结</h2>
<ul>
<li>文件：<code>with open</code> + 逐行读取处理大文件；CSV 用 <code>csv</code> 模块或 Pandas</li>
<li>正则：<code>re.search</code> / <code>findall</code> / <code>sub</code>，先写小模式再组合</li>
<li>NumPy：向量化运算 + 广播，避免 Python 循环</li>
<li>Pandas：<code>read_csv</code> → 筛选/分组/缺失值处理 → 导出，是表格数据的标准工作流</li>
</ul>
<p>至此 Python 教程三部曲完成：入门 → 进阶 → 数据处理。</p>`,E_=`<h1>R Language Primer: Data Structures, Vectorization, and Functions</h1>
<p>R is the go-to language for statistical computing and data visualization. This article explains the core of R: five data structures, vectorized operations, control flow, and functions, laying the groundwork for subsequent tidyverse and ggplot2.</p>
<h2>1. Environment Setup</h2>
<h3>1.1 Installation and IDE</h3>
<ul>
<li>Install <a href="https://cran.r-project.org/">R</a></li>
<li>Recommended IDE: <a href="https://posit.co/download/rstudio-desktop/">RStudio</a> (free), or VS Code + R extension</li>
</ul>
<h3>1.2 Basic Operations</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Assignment: &#x3C;- is the classic R assignment operator (= also works)</span>
x <span class="hljs-operator">&#x3C;-</span> 10
x

<span class="hljs-comment"># View help</span>
<span class="hljs-operator">?</span>mean
help<span class="hljs-punctuation">(</span><span class="hljs-string">"lm"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>2. Basic Data Structures</h2>
<p>R's core data structures can be classified by dimension and homogeneity:</p>
<table>
<thead>
<tr>
<th>Structure</th>
<th>Dimension</th>
<th>Element Type</th>
</tr>
</thead>
<tbody>
<tr>
<td>vector</td>
<td>1</td>
<td>Homogeneous</td>
</tr>
<tr>
<td>factor</td>
<td>1</td>
<td>Categorical</td>
</tr>
<tr>
<td>matrix</td>
<td>2</td>
<td>Homogeneous</td>
</tr>
<tr>
<td>data.frame</td>
<td>2</td>
<td>Heterogeneous allowed</td>
</tr>
<tr>
<td>list</td>
<td>N</td>
<td>Heterogeneous allowed</td>
</tr>
</tbody>
</table>
<h3>2.1 Vectors</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Create vectors</span>
a <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>            <span class="hljs-comment"># combine</span>
b <span class="hljs-operator">&#x3C;-</span> 1<span class="hljs-operator">:</span><span class="hljs-number">10</span>                        <span class="hljs-comment"># sequence</span>
<span class="hljs-built_in">c</span> <span class="hljs-operator">&#x3C;-</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">)</span>         <span class="hljs-comment"># step sequence</span>
d <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-string">"A"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>                 <span class="hljs-comment"># repeat</span>
e <span class="hljs-operator">&#x3C;-</span> sample<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">100</span><span class="hljs-punctuation">,</span> <span class="hljs-number">10</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># random sample</span>

<span class="hljs-comment"># Types</span>
typeof<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># "double"</span>
<span class="hljs-built_in">is.numeric</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>  <span class="hljs-comment"># TRUE</span>

<span class="hljs-comment"># Indexing (R starts at 1!)</span>
a<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span>           <span class="hljs-comment"># 1</span>
a<span class="hljs-punctuation">[</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>     <span class="hljs-comment"># elements 1 and 3</span>
a<span class="hljs-punctuation">[</span><span class="hljs-operator">-</span><span class="hljs-number">2</span><span class="hljs-punctuation">]</span>          <span class="hljs-comment"># remove element 2</span>
a<span class="hljs-punctuation">[</span>a <span class="hljs-operator">></span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>       <span class="hljs-comment"># conditional filtering</span>
<span class="hljs-built_in">names</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span> <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x2"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x3"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x4"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x5"</span><span class="hljs-punctuation">)</span>
a<span class="hljs-punctuation">[</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># access by name</span>
</code></pre>
<blockquote>
<p><strong>R indexes from 1</strong>, which is the biggest habit difference from Python.</p>
</blockquote>
<h3>2.2 Factors</h3>
<pre><code class="hljs language-r">treatment <span class="hljs-operator">&#x3C;-</span> factor<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
                    levels <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
levels<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>          <span class="hljs-comment"># "control" "drug"</span>
table<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># frequency counts</span>
</code></pre>
<h3>2.3 Matrices</h3>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> ncol <span class="hljs-operator">=</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment">#      [,1] [,2] [,3] [,4]</span>
<span class="hljs-comment"># [1,]    1    4    7   10</span>
<span class="hljs-comment"># [2,]    2    5    8   11</span>
<span class="hljs-comment"># [3,]    3    6    9   12</span>

m<span class="hljs-punctuation">[</span><span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># row 2, column 3: 8</span>
m<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>         <span class="hljs-comment"># row 1</span>
m<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">]</span>         <span class="hljs-comment"># column 2</span>
t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># transpose</span>
m <span class="hljs-operator">%*%</span> t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>     <span class="hljs-comment"># matrix multiplication</span>
</code></pre>
<h3>2.4 Data Frames</h3>
<p>Data frames are the core structure for tabular data:</p>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># View structure</span>
str<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>
summary<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Extract columns</span>
df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span>          <span class="hljs-comment"># $ operator</span>
df<span class="hljs-punctuation">[[</span><span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>
df<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span>

<span class="hljs-comment"># Extract rows</span>
df<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>                <span class="hljs-comment"># first row</span>
df<span class="hljs-punctuation">[</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>

<span class="hljs-comment"># Add column</span>
df<span class="hljs-operator">$</span>log2fc <span class="hljs-operator">&#x3C;-</span> log2<span class="hljs-punctuation">(</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.5 Lists</h3>
<p>Lists can hold elements of any type:</p>
<pre><code class="hljs language-r">result <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>
  name <span class="hljs-operator">=</span> <span class="hljs-string">"Analysis Result"</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> <span class="hljs-number">0.003</span><span class="hljs-punctuation">,</span>
  coef <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1.2</span><span class="hljs-punctuation">,</span> <span class="hljs-operator">-</span><span class="hljs-number">0.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  model <span class="hljs-operator">=</span> lm<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">~</span> chromosome<span class="hljs-punctuation">,</span> data <span class="hljs-operator">=</span> df<span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

result<span class="hljs-operator">$</span>pvalue
result<span class="hljs-punctuation">[[</span><span class="hljs-number">3</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>            <span class="hljs-comment"># numeric vector</span>
</code></pre>
<h2>3. Vectorized Operations</h2>
<p><strong>R's core idea: operate on entire vectors, not element-by-element loops.</strong></p>
<pre><code class="hljs language-r">x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>

x <span class="hljs-operator">*</span> <span class="hljs-number">2</span>              <span class="hljs-comment"># element-wise multiplication</span>
x <span class="hljs-operator">+</span> <span class="hljs-number">10</span>
<span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
log2<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># vector with vector</span>
y <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">+</span> y              <span class="hljs-comment"># add corresponding positions</span>

<span class="hljs-comment"># comparison operations return logical vectors</span>
x <span class="hljs-operator">></span> <span class="hljs-number">3</span>              <span class="hljs-comment"># FALSE FALSE FALSE  TRUE  TRUE</span>

<span class="hljs-comment"># statistical functions</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; median<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; <span class="hljs-built_in">range</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. Control Flow</h2>
<h3>4.1 Conditionals</h3>
<pre><code class="hljs language-r">score <span class="hljs-operator">&#x3C;-</span> 85

<span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">90</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"A"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">80</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"B"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"C"</span>
<span class="hljs-punctuation">}</span>
print<span class="hljs-punctuation">(</span>grade<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># vectorized conditional: ifelse</span>
values <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45</span><span class="hljs-punctuation">,</span> <span class="hljs-number">75</span><span class="hljs-punctuation">)</span>
result <span class="hljs-operator">&#x3C;-</span> ifelse<span class="hljs-punctuation">(</span>values <span class="hljs-operator">>=</span> <span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"pass"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"fail"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>result<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># "pass" "pass" "fail" "pass"</span>
</code></pre>
<h3>4.2 Loops</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># for loop</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>i <span class="hljs-keyword">in</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># iterate over vector elements</span>
genes <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>g <span class="hljs-keyword">in</span> genes<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>paste<span class="hljs-punctuation">(</span><span class="hljs-string">"Analyzing"</span><span class="hljs-punctuation">,</span> g<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># while loop</span>
n <span class="hljs-operator">&#x3C;-</span> 0
<span class="hljs-keyword">while</span> <span class="hljs-punctuation">(</span>n <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  n <span class="hljs-operator">&#x3C;-</span> n <span class="hljs-operator">+</span> <span class="hljs-number">1</span>
<span class="hljs-punctuation">}</span>
</code></pre>
<h3>4.3 Avoid Loops When Possible: The apply Family</h3>
<p>R commonly uses the <code>apply</code> family instead of explicit loops:</p>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># row-wise means</span>
apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">)</span>       <span class="hljs-comment"># column-wise sums</span>

lapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># apply function to each element of a list, returns a list</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># simplified version, returns a vector</span>
</code></pre>
<h2>5. Functions</h2>
<h3>5.1 Defining Functions</h3>
<pre><code class="hljs language-r">gc_content <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  seq <span class="hljs-operator">&#x3C;-</span> toupper<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>strsplit<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">,</span> <span class="hljs-string">""</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">[[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span> <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"G"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"C"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">/</span> nchar<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-operator">*</span> <span class="hljs-number">100</span>
<span class="hljs-punctuation">}</span>

gc_content<span class="hljs-punctuation">(</span><span class="hljs-string">"ATGCCGA"</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment"># [1] 57.14286</span>
</code></pre>
<h3>5.2 Default Arguments and Return Values</h3>
<pre><code class="hljs language-r">normalize <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> <span class="hljs-punctuation">(</span><span class="hljs-built_in">max</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span> <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
    stop<span class="hljs-punctuation">(</span><span class="hljs-string">"Unknown method: "</span><span class="hljs-punctuation">,</span> method<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span>
<span class="hljs-punctuation">}</span>

normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>5.3 Anonymous Functions and Functional Programming</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># anonymous functions</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>i<span class="hljs-punctuation">)</span> i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># purrr-style map (requires tidyverse, detailed in the next article)</span>
<span class="hljs-comment"># purrr::map_dbl(1:5, ~ .x ^ 2)</span>
</code></pre>
<h2>6. Practical Tips</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># pipe (native pipe |>/ in R 4.1+ or magrittr %>%)</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># handle missing values</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">)</span>
<span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
x<span class="hljs-punctuation">[</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>
na.omit<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># merging data frames</span>
df1 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> value <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"a"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"b"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"c"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
df2 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-operator">:</span><span class="hljs-number">4</span><span class="hljs-punctuation">,</span> score <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">85</span><span class="hljs-punctuation">,</span> <span class="hljs-number">88</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
merged <span class="hljs-operator">&#x3C;-</span> merge<span class="hljs-punctuation">(</span>df1<span class="hljs-punctuation">,</span> df2<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"id"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># similar to SQL join</span>
</code></pre>
<h2>7. Summary</h2>
<ul>
<li>Five core structures: vectors / factors / matrices / data frames / lists</li>
<li><strong>R indexing starts at 1</strong>; <code>$</code> for columns; <code>[ ]</code> for rows</li>
<li>Vectorized operations + <code>ifelse</code> + the <code>apply</code> family replace explicit loops</li>
<li>Functions: default arguments, <code>stop()</code> for errors, anonymous functions</li>
</ul>
<p>The next article will introduce tidyverse: using <code>dplyr</code> for elegant data manipulation.</p>`,T_=`<h1>R 语言入门：数据结构、向量化与函数</h1>
<p>R 是统计计算与数据可视化的首选语言。本文讲解 R 的核心：五种数据结构、向量化运算思维、控制流与函数，为后续 tidyverse 与 ggplot2 打下基础。</p>
<h2>1. 环境准备</h2>
<h3>1.1 安装与 IDE</h3>
<ul>
<li>安装 <a href="https://cran.r-project.org/">R</a></li>
<li>推荐 IDE：<a href="https://posit.co/download/rstudio-desktop/">RStudio</a>（免费），或 VS Code + R 扩展</li>
</ul>
<h3>1.2 基本操作</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 赋值：&#x3C;- 是 R 的经典赋值符（= 也可用）</span>
x <span class="hljs-operator">&#x3C;-</span> 10
x

<span class="hljs-comment"># 查看帮助</span>
<span class="hljs-operator">?</span>mean
help<span class="hljs-punctuation">(</span><span class="hljs-string">"lm"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>2. 基础数据结构</h2>
<p>R 的核心数据结构按维度与同质性划分：</p>
<table>
<thead>
<tr>
<th>结构</th>
<th>维度</th>
<th>元素类型</th>
</tr>
</thead>
<tbody>
<tr>
<td>vector（向量）</td>
<td>1</td>
<td>同质</td>
</tr>
<tr>
<td>factor（因子）</td>
<td>1</td>
<td>分类</td>
</tr>
<tr>
<td>matrix（矩阵）</td>
<td>2</td>
<td>同质</td>
</tr>
<tr>
<td>data.frame（数据框）</td>
<td>2</td>
<td>可异质</td>
</tr>
<tr>
<td>list（列表）</td>
<td>N</td>
<td>可异质</td>
</tr>
</tbody>
</table>
<h3>2.1 向量（vector）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 创建向量</span>
a <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>            <span class="hljs-comment"># combine</span>
b <span class="hljs-operator">&#x3C;-</span> 1<span class="hljs-operator">:</span><span class="hljs-number">10</span>                        <span class="hljs-comment"># 序列</span>
<span class="hljs-built_in">c</span> <span class="hljs-operator">&#x3C;-</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">)</span>         <span class="hljs-comment"># 步长序列</span>
d <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-string">"A"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>                 <span class="hljs-comment"># 重复</span>
e <span class="hljs-operator">&#x3C;-</span> sample<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">100</span><span class="hljs-punctuation">,</span> <span class="hljs-number">10</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 随机抽样</span>

<span class="hljs-comment"># 类型</span>
typeof<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># "double"</span>
<span class="hljs-built_in">is.numeric</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>  <span class="hljs-comment"># TRUE</span>

<span class="hljs-comment"># 索引（R 从 1 开始！）</span>
a<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span>           <span class="hljs-comment"># 1</span>
a<span class="hljs-punctuation">[</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>     <span class="hljs-comment"># 第 1、3 个元素</span>
a<span class="hljs-punctuation">[</span><span class="hljs-operator">-</span><span class="hljs-number">2</span><span class="hljs-punctuation">]</span>          <span class="hljs-comment"># 删除第 2 个</span>
a<span class="hljs-punctuation">[</span>a <span class="hljs-operator">></span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>       <span class="hljs-comment"># 条件筛选</span>
<span class="hljs-built_in">names</span><span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span> <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x2"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x3"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x4"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"x5"</span><span class="hljs-punctuation">)</span>
a<span class="hljs-punctuation">[</span><span class="hljs-string">"x1"</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># 按名字取</span>
</code></pre>
<blockquote>
<p><strong>R 的索引从 1 开始</strong>，这是与 Python 最大的习惯差异。</p>
</blockquote>
<h3>2.2 因子（factor）</h3>
<pre><code class="hljs language-r">treatment <span class="hljs-operator">&#x3C;-</span> factor<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
                    levels <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"drug"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
levels<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>          <span class="hljs-comment"># "control" "drug"</span>
table<span class="hljs-punctuation">(</span>treatment<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 频数统计</span>
</code></pre>
<h3>2.3 矩阵（matrix）</h3>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> ncol <span class="hljs-operator">=</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment">#      [,1] [,2] [,3] [,4]</span>
<span class="hljs-comment"># [1,]    1    4    7   10</span>
<span class="hljs-comment"># [2,]    2    5    8   11</span>
<span class="hljs-comment"># [3,]    3    6    9   12</span>

m<span class="hljs-punctuation">[</span><span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">]</span>        <span class="hljs-comment"># 第 2 行第 3 列：8</span>
m<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>         <span class="hljs-comment"># 第 1 行</span>
m<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">]</span>         <span class="hljs-comment"># 第 2 列</span>
t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 转置</span>
m <span class="hljs-operator">%*%</span> t<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">)</span>     <span class="hljs-comment"># 矩阵乘法</span>
</code></pre>
<h3>2.4 数据框（data.frame）</h3>
<p>数据框是表格数据的核心结构：</p>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 查看结构</span>
str<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>
summary<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 取列</span>
df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span>          <span class="hljs-comment"># $ 运算符</span>
df<span class="hljs-punctuation">[[</span><span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>
df<span class="hljs-punctuation">[</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">]</span>

<span class="hljs-comment"># 取行</span>
df<span class="hljs-punctuation">[</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>                <span class="hljs-comment"># 第一行</span>
df<span class="hljs-punctuation">[</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span><span class="hljs-punctuation">,</span> <span class="hljs-punctuation">]</span>

<span class="hljs-comment"># 增加列</span>
df<span class="hljs-operator">$</span>log2fc <span class="hljs-operator">&#x3C;-</span> log2<span class="hljs-punctuation">(</span>df<span class="hljs-operator">$</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.5 列表（list）</h3>
<p>列表可容纳任意类型：</p>
<pre><code class="hljs language-r">result <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>
  name <span class="hljs-operator">=</span> <span class="hljs-string">"分析结果"</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> <span class="hljs-number">0.003</span><span class="hljs-punctuation">,</span>
  coef <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1.2</span><span class="hljs-punctuation">,</span> <span class="hljs-operator">-</span><span class="hljs-number">0.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  model <span class="hljs-operator">=</span> lm<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">~</span> chromosome<span class="hljs-punctuation">,</span> data <span class="hljs-operator">=</span> df<span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

result<span class="hljs-operator">$</span>pvalue
result<span class="hljs-punctuation">[[</span><span class="hljs-number">3</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span>            <span class="hljs-comment"># 数值向量</span>
</code></pre>
<h2>3. 向量化运算</h2>
<p><strong>R 的核心思想：对整个向量运算，而非逐元素循环。</strong></p>
<pre><code class="hljs language-r">x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>

x <span class="hljs-operator">*</span> <span class="hljs-number">2</span>              <span class="hljs-comment"># 逐元素乘</span>
x <span class="hljs-operator">+</span> <span class="hljs-number">10</span>
<span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
log2<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 向量与向量</span>
y <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">+</span> y              <span class="hljs-comment"># 对应位置相加</span>

<span class="hljs-comment"># 比较运算返回逻辑向量</span>
x <span class="hljs-operator">></span> <span class="hljs-number">3</span>              <span class="hljs-comment"># FALSE FALSE FALSE  TRUE  TRUE</span>

<span class="hljs-comment"># 统计函数</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; median<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>; <span class="hljs-built_in">range</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. 控制流</h2>
<h3>4.1 条件</h3>
<pre><code class="hljs language-r">score <span class="hljs-operator">&#x3C;-</span> 85

<span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">90</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"A"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>score <span class="hljs-operator">>=</span> <span class="hljs-number">80</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"B"</span>
<span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
  grade <span class="hljs-operator">&#x3C;-</span> <span class="hljs-string">"C"</span>
<span class="hljs-punctuation">}</span>
print<span class="hljs-punctuation">(</span>grade<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 向量化条件：ifelse</span>
values <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45</span><span class="hljs-punctuation">,</span> <span class="hljs-number">75</span><span class="hljs-punctuation">)</span>
result <span class="hljs-operator">&#x3C;-</span> ifelse<span class="hljs-punctuation">(</span>values <span class="hljs-operator">>=</span> <span class="hljs-number">60</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"pass"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"fail"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>result<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># "pass" "pass" "fail" "pass"</span>
</code></pre>
<h3>4.2 循环</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># for 循环</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>i <span class="hljs-keyword">in</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># 遍历向量元素</span>
genes <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">)</span>
<span class="hljs-keyword">for</span> <span class="hljs-punctuation">(</span>g <span class="hljs-keyword">in</span> genes<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  print<span class="hljs-punctuation">(</span>paste<span class="hljs-punctuation">(</span><span class="hljs-string">"分析"</span><span class="hljs-punctuation">,</span> g<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">}</span>

<span class="hljs-comment"># while 循环</span>
n <span class="hljs-operator">&#x3C;-</span> 0
<span class="hljs-keyword">while</span> <span class="hljs-punctuation">(</span>n <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  n <span class="hljs-operator">&#x3C;-</span> n <span class="hljs-operator">+</span> <span class="hljs-number">1</span>
<span class="hljs-punctuation">}</span>
</code></pre>
<h3>4.3 尽量避免循环：apply 家族</h3>
<p>R 中常用 <code>apply</code> 家族替代显式循环：</p>
<pre><code class="hljs language-r">m <span class="hljs-operator">&#x3C;-</span> matrix<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">12</span><span class="hljs-punctuation">,</span> nrow <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>      <span class="hljs-comment"># 对每行求均值</span>
apply<span class="hljs-punctuation">(</span>m<span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">)</span>       <span class="hljs-comment"># 对每列求和</span>

lapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 对列表每个元素操作，返回列表</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-operator">:</span><span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> mean<span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 简化版，返回向量</span>
</code></pre>
<h2>5. 函数</h2>
<h3>5.1 定义函数</h3>
<pre><code class="hljs language-r">gc_content <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  seq <span class="hljs-operator">&#x3C;-</span> toupper<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span>strsplit<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">,</span> <span class="hljs-string">""</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">[[</span><span class="hljs-number">1</span><span class="hljs-punctuation">]</span><span class="hljs-punctuation">]</span> <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"G"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"C"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  gc <span class="hljs-operator">/</span> nchar<span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">)</span> <span class="hljs-operator">*</span> <span class="hljs-number">100</span>
<span class="hljs-punctuation">}</span>

gc_content<span class="hljs-punctuation">(</span><span class="hljs-string">"ATGCCGA"</span><span class="hljs-punctuation">)</span>
<span class="hljs-comment"># [1] 57.14286</span>
</code></pre>
<h3>5.2 默认参数与返回值</h3>
<pre><code class="hljs language-r">normalize <span class="hljs-operator">&#x3C;-</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
  <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"zscore"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> mean<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> sd<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-keyword">if</span> <span class="hljs-punctuation">(</span>method <span class="hljs-operator">==</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span> <span class="hljs-punctuation">{</span>
    <span class="hljs-punctuation">(</span>x <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">/</span> <span class="hljs-punctuation">(</span><span class="hljs-built_in">max</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span> <span class="hljs-operator">-</span> <span class="hljs-built_in">min</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span> <span class="hljs-keyword">else</span> <span class="hljs-punctuation">{</span>
    stop<span class="hljs-punctuation">(</span><span class="hljs-string">"未知方法: "</span><span class="hljs-punctuation">,</span> method<span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">}</span>
<span class="hljs-punctuation">}</span>

normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
normalize<span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> method <span class="hljs-operator">=</span> <span class="hljs-string">"minmax"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>5.3 匿名函数与函数式编程</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 匿名函数</span>
sapply<span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-keyword">function</span><span class="hljs-punctuation">(</span>i<span class="hljs-punctuation">)</span> i <span class="hljs-operator">^</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># purrr 风格的 map（需要 tidyverse，下一篇详述）</span>
<span class="hljs-comment"># purrr::map_dbl(1:5, ~ .x ^ 2)</span>
</code></pre>
<h2>6. 实用技巧</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 管道（R 4.1+ 原生管道 |>/ 或 magrittr %>%）</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">)</span>
x <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 缺失值处理</span>
x <span class="hljs-operator">&#x3C;-</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">)</span>
<span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>
<span class="hljs-built_in">sum</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
x<span class="hljs-punctuation">[</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span><span class="hljs-punctuation">]</span>
na.omit<span class="hljs-punctuation">(</span>x<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 合并数据框</span>
df1 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> value <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"a"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"b"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"c"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
df2 <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>id <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-operator">:</span><span class="hljs-number">4</span><span class="hljs-punctuation">,</span> score <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">90</span><span class="hljs-punctuation">,</span> <span class="hljs-number">85</span><span class="hljs-punctuation">,</span> <span class="hljs-number">88</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
merged <span class="hljs-operator">&#x3C;-</span> merge<span class="hljs-punctuation">(</span>df1<span class="hljs-punctuation">,</span> df2<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"id"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 类似 SQL join</span>
</code></pre>
<h2>7. 小结</h2>
<ul>
<li>五种核心结构：向量 / 因子 / 矩阵 / 数据框 / 列表</li>
<li><strong>R 索引从 1 开始</strong>；<code>$</code> 取列；<code>[ ]</code> 取行</li>
<li>向量化运算 + <code>ifelse</code> + <code>apply</code> 家族替代显式循环</li>
<li>函数：默认参数、<code>stop()</code> 报错、匿名函数</li>
</ul>
<p>下一篇将介绍 tidyverse：用 <code>dplyr</code> 优雅地完成数据操作。</p>`,A_=`<h1>R ggplot2 Data Visualization: Grammar of Graphics and Common Charts</h1>
<p>ggplot2 is R's most powerful plotting package, built on the <strong>Grammar of Graphics</strong>: data, mapping, geometric objects, statistical transformations, coordinates, facets, and themes, stacked layer by layer into a figure.</p>
<h2>1. Basic Principles</h2>
<pre><code class="hljs language-r">library<span class="hljs-punctuation">(</span>ggplot2<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Skeleton of a plot</span>
ggplot<span class="hljs-punctuation">(</span>data <span class="hljs-operator">=</span> dataset<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> variable<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> variable<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>          <span class="hljs-comment"># Geometric layer (point/line/bar/boxplot...)</span>
  labs<span class="hljs-punctuation">(</span>...<span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>           <span class="hljs-comment"># Labels</span>
  theme_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># Theme</span>
</code></pre>
<ul>
<li><code>data</code>: data frame</li>
<li><code>aes()</code>: aesthetic mapping (x, y, color, fill, size, shape)</li>
<li><code>geom_*</code>: geometric objects, one per layer</li>
</ul>
<h2>2. Quick Start with Built-in Datasets</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># mpg: car fuel economy data; iris: iris flower data</span>
head<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">)</span>
head<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.1 Scatter Plot</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Color by category + smooth trend line</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.7</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_smooth<span class="hljs-punctuation">(</span>method <span class="hljs-operator">=</span> <span class="hljs-string">"lm"</span><span class="hljs-punctuation">,</span> se <span class="hljs-operator">=</span> <span class="hljs-literal">FALSE</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"black"</span><span class="hljs-punctuation">,</span> linewidth <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.2 Boxplot</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_boxplot<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_jitter<span class="hljs-punctuation">(</span>width <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.5</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># Overlay points to show distribution</span>
</code></pre>
<h3>2.3 Histogram</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Petal.Length<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_histogram<span class="hljs-punctuation">(</span>bins <span class="hljs-operator">=</span> <span class="hljs-number">30</span><span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> <span class="hljs-string">"steelblue"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"white"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.4 Bar Chart</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Count bar chart</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_bar<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Value bar chart (requires aggregation first)</span>
library<span class="hljs-punctuation">(</span>dplyr<span class="hljs-punctuation">)</span>
avg_hwy <span class="hljs-operator">&#x3C;-</span> mpg <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_hwy <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>avg_hwy<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> reorder<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">,</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_col<span class="hljs-punctuation">(</span>fill <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  coord_flip<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>          <span class="hljs-comment"># Horizontal, easier to read long category names</span>
</code></pre>
<h3>2.5 Line Chart</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Time series example</span>
df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  time <span class="hljs-operator">=</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  value <span class="hljs-operator">=</span> <span class="hljs-built_in">sin</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">11</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> time<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> value<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_line<span class="hljs-punctuation">(</span>linewidth <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>3. Faceting (facet)</h2>
<p>Split into multiple subplots by variables:</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># Facet by drv: one row</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_wrap<span class="hljs-punctuation">(</span><span class="hljs-operator">~</span>drv<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Two-variable faceting: grid</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_grid<span class="hljs-punctuation">(</span>cyl <span class="hljs-operator">~</span> drv<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. Color Themes</h2>
<h3>4.1 Scale Control</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Continuous color scale</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> year<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_gradient<span class="hljs-punctuation">(</span>low <span class="hljs-operator">=</span> <span class="hljs-string">"blue"</span><span class="hljs-punctuation">,</span> high <span class="hljs-operator">=</span> <span class="hljs-string">"red"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Discrete color palette</span>
ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_brewer<span class="hljs-punctuation">(</span>palette <span class="hljs-operator">=</span> <span class="hljs-string">"Set1"</span><span class="hljs-punctuation">)</span>     <span class="hljs-comment"># Built-in RColorBrewer palette</span>

<span class="hljs-comment"># Manually specified</span>
scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#10b981"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#f59e0b"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 Theme Customization</h3>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Classic white background theme</span>
p <span class="hljs-operator">+</span> theme_bw<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Minimal theme</span>
p <span class="hljs-operator">+</span> theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Custom details</span>
p <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span>base_size <span class="hljs-operator">=</span> <span class="hljs-number">14</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme<span class="hljs-punctuation">(</span>
    legend.position <span class="hljs-operator">=</span> <span class="hljs-string">"top"</span><span class="hljs-punctuation">,</span>                 <span class="hljs-comment"># Legend position</span>
    panel.grid.minor <span class="hljs-operator">=</span> element_blank<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># Remove minor grid lines</span>
    plot.title <span class="hljs-operator">=</span> element_text<span class="hljs-punctuation">(</span>face <span class="hljs-operator">=</span> <span class="hljs-string">"bold"</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">16</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. Labels and Export</h2>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>
    title <span class="hljs-operator">=</span> <span class="hljs-string">"Iris Measurement Data"</span><span class="hljs-punctuation">,</span>
    subtitle <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal size distribution"</span><span class="hljs-punctuation">,</span>
    x <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal Length (cm)"</span><span class="hljs-punctuation">,</span>
    y <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal Width (cm)"</span><span class="hljs-punctuation">,</span>
    color <span class="hljs-operator">=</span> <span class="hljs-string">"Species"</span>
  <span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Save: ggsave automatically recognizes the extension</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.png"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">,</span> dpi <span class="hljs-operator">=</span> <span class="hljs-number">300</span><span class="hljs-punctuation">)</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.pdf"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span>
</code></pre>
<blockquote>
<p>For journal submission, export <strong>PDF vector graphics</strong> (300 dpi PNG for web).</p>
</blockquote>
<h2>6. Bioinformatics Practice: Volcano Plot</h2>
<p>Standard visualization for differential expression analysis:</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># Simulate differential expression results</span>
set.seed<span class="hljs-punctuation">(</span><span class="hljs-number">42</span><span class="hljs-punctuation">)</span>
de <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> paste0<span class="hljs-punctuation">(</span><span class="hljs-string">"GENE"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">500</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  log2fc <span class="hljs-operator">=</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1.8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> runif<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Mark significant genes</span>
de <span class="hljs-operator">&#x3C;-</span> de <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>significance <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">abs</span><span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span> <span class="hljs-operator">></span> <span class="hljs-number">1</span> <span class="hljs-operator">&#x26;</span> pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"up/down"</span><span class="hljs-punctuation">,</span>
    pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"significant"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span> <span class="hljs-operator">~</span> <span class="hljs-string">"ns"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>de<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> log2fc<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span>pvalue<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> significance<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">1.5</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"up/down"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#ef4444"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"significant"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"ns"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#9ca3af"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_vline<span class="hljs-punctuation">(</span>xintercept <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-operator">-</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_hline<span class="hljs-punctuation">(</span>yintercept <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span><span class="hljs-number">0.05</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>title <span class="hljs-operator">=</span> <span class="hljs-string">"Differential Expression Volcano Plot"</span><span class="hljs-punctuation">,</span> x <span class="hljs-operator">=</span> <span class="hljs-string">"log2 Fold Change"</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-string">"-log10(p-value)"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>7. Summary</h2>
<ul>
<li>Grammar of Graphics: <code>ggplot(data, aes()) + geom_xxx() + labs() + theme_xxx()</code></li>
<li>Common geometries: <code>geom_point</code> / <code>geom_boxplot</code> / <code>geom_histogram</code> / <code>geom_col</code> / <code>geom_line</code></li>
<li>Group by <code>aes(color/fill)</code>, split plots by <code>facet_wrap</code> / <code>facet_grid</code></li>
<li><code>theme_minimal</code> + <code>labs</code> + <code>ggsave(dpi=300)</code> delivers publication-ready output</li>
</ul>
<p>With this, the R tutorial trilogy is complete: Introduction → tidyverse → ggplot2.</p>`,R_=`<h1>R ggplot2 数据可视化：图层语法与常用图表</h1>
<p>ggplot2 是 R 最强大的绘图包，基于<strong>图层语法</strong>（Grammar of Graphics）：数据、映射、几何对象、统计变换、坐标、分面、主题，逐层叠加成图。</p>
<h2>1. 基本原理</h2>
<pre><code class="hljs language-r">library<span class="hljs-punctuation">(</span>ggplot2<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 一张图的骨架</span>
ggplot<span class="hljs-punctuation">(</span>data <span class="hljs-operator">=</span> 数据集<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> 变量<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> 变量<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>          <span class="hljs-comment"># 几何图层（点/线/柱/箱线...）</span>
  labs<span class="hljs-punctuation">(</span>...<span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>           <span class="hljs-comment"># 标签</span>
  theme_xxx<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 主题</span>
</code></pre>
<ul>
<li><code>data</code>：数据框</li>
<li><code>aes()</code>：美学映射（x、y、color、fill、size、shape）</li>
<li><code>geom_*</code>：几何对象，每层一个</li>
</ul>
<h2>2. 内置数据集快速上手</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># mpg：汽车油耗数据；iris：鸢尾花数据</span>
head<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">)</span>
head<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.1 散点图（scatter）</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 按类别着色 + 平滑趋势线</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.7</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_smooth<span class="hljs-punctuation">(</span>method <span class="hljs-operator">=</span> <span class="hljs-string">"lm"</span><span class="hljs-punctuation">,</span> se <span class="hljs-operator">=</span> <span class="hljs-literal">FALSE</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"black"</span><span class="hljs-punctuation">,</span> linewidth <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.2 箱线图（boxplot）</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_boxplot<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_jitter<span class="hljs-punctuation">(</span>width <span class="hljs-operator">=</span> <span class="hljs-number">0.2</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.5</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 叠加散点展示分布</span>
</code></pre>
<h3>2.3 直方图（histogram）</h3>
<pre><code class="hljs language-r">ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Petal.Length<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_histogram<span class="hljs-punctuation">(</span>bins <span class="hljs-operator">=</span> <span class="hljs-number">30</span><span class="hljs-punctuation">,</span> fill <span class="hljs-operator">=</span> <span class="hljs-string">"steelblue"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"white"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>2.4 柱状图（bar）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 计数柱状图</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> <span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_bar<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 数值柱状图（需要先聚合）</span>
library<span class="hljs-punctuation">(</span>dplyr<span class="hljs-punctuation">)</span>
avg_hwy <span class="hljs-operator">&#x3C;-</span> mpg <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_hwy <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>avg_hwy<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> reorder<span class="hljs-punctuation">(</span><span class="hljs-built_in">class</span><span class="hljs-punctuation">,</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> mean_hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_col<span class="hljs-punctuation">(</span>fill <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  coord_flip<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>          <span class="hljs-comment"># 横置，长类别名更易读</span>
</code></pre>
<h3>2.5 折线图（line）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 时间序列示例</span>
df <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  time <span class="hljs-operator">=</span> seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  value <span class="hljs-operator">=</span> <span class="hljs-built_in">sin</span><span class="hljs-punctuation">(</span>seq<span class="hljs-punctuation">(</span><span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">20</span><span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">11</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>df<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> time<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> value<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_line<span class="hljs-punctuation">(</span>linewidth <span class="hljs-operator">=</span> <span class="hljs-number">1</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>3. 分面（facet）</h2>
<p>按变量拆分成多个子图：</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># 按 drv 分面：一排</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_wrap<span class="hljs-punctuation">(</span><span class="hljs-operator">~</span>drv<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 双变量分面：网格</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  facet_grid<span class="hljs-punctuation">(</span>cyl <span class="hljs-operator">~</span> drv<span class="hljs-punctuation">)</span>
</code></pre>
<h2>4. 颜色主题</h2>
<h3>4.1 标度控制（scale）</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 连续色阶</span>
ggplot<span class="hljs-punctuation">(</span>mpg<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> displ<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> hwy<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> year<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_gradient<span class="hljs-punctuation">(</span>low <span class="hljs-operator">=</span> <span class="hljs-string">"blue"</span><span class="hljs-punctuation">,</span> high <span class="hljs-operator">=</span> <span class="hljs-string">"red"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 离散色板</span>
ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_brewer<span class="hljs-punctuation">(</span>palette <span class="hljs-operator">=</span> <span class="hljs-string">"Set1"</span><span class="hljs-punctuation">)</span>     <span class="hljs-comment"># 内置 RColorBrewer 色板</span>

<span class="hljs-comment"># 手动指定</span>
scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#10b981"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"#f59e0b"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 主题美化</h3>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 经典白底主题</span>
p <span class="hljs-operator">+</span> theme_bw<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 极简主题</span>
p <span class="hljs-operator">+</span> theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 自定义细节</span>
p <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span>base_size <span class="hljs-operator">=</span> <span class="hljs-number">14</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme<span class="hljs-punctuation">(</span>
    legend.position <span class="hljs-operator">=</span> <span class="hljs-string">"top"</span><span class="hljs-punctuation">,</span>                 <span class="hljs-comment"># 图例位置</span>
    panel.grid.minor <span class="hljs-operator">=</span> element_blank<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># 去掉次要网格线</span>
    plot.title <span class="hljs-operator">=</span> element_text<span class="hljs-punctuation">(</span>face <span class="hljs-operator">=</span> <span class="hljs-string">"bold"</span><span class="hljs-punctuation">,</span> size <span class="hljs-operator">=</span> <span class="hljs-number">16</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. 标签与导出</h2>
<pre><code class="hljs language-r">p <span class="hljs-operator">&#x3C;-</span> ggplot<span class="hljs-punctuation">(</span>iris<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> Sepal.Length<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> Sepal.Width<span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> Species<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">2.5</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>
    title <span class="hljs-operator">=</span> <span class="hljs-string">"鸢尾花测量数据"</span><span class="hljs-punctuation">,</span>
    subtitle <span class="hljs-operator">=</span> <span class="hljs-string">"Sepal 尺寸分布"</span><span class="hljs-punctuation">,</span>
    x <span class="hljs-operator">=</span> <span class="hljs-string">"花萼长度 (cm)"</span><span class="hljs-punctuation">,</span>
    y <span class="hljs-operator">=</span> <span class="hljs-string">"花萼宽度 (cm)"</span><span class="hljs-punctuation">,</span>
    color <span class="hljs-operator">=</span> <span class="hljs-string">"物种"</span>
  <span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 保存：ggsave 自动识别扩展名</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.png"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">,</span> dpi <span class="hljs-operator">=</span> <span class="hljs-number">300</span><span class="hljs-punctuation">)</span>
ggsave<span class="hljs-punctuation">(</span><span class="hljs-string">"iris_scatter.pdf"</span><span class="hljs-punctuation">,</span> p<span class="hljs-punctuation">,</span> width <span class="hljs-operator">=</span> <span class="hljs-number">8</span><span class="hljs-punctuation">,</span> height <span class="hljs-operator">=</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span>
</code></pre>
<blockquote>
<p>论文投稿建议导出 <strong>PDF 矢量图</strong>（300 dpi 的 PNG 用于网页）。</p>
</blockquote>
<h2>6. 生物信息学实战：火山图</h2>
<p>差异表达分析的标准可视化：</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># 模拟差异表达结果</span>
set.seed<span class="hljs-punctuation">(</span><span class="hljs-number">42</span><span class="hljs-punctuation">)</span>
de <span class="hljs-operator">&#x3C;-</span> data.frame<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> paste0<span class="hljs-punctuation">(</span><span class="hljs-string">"GENE"</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">500</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  log2fc <span class="hljs-operator">=</span> rnorm<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1.8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  pvalue <span class="hljs-operator">=</span> runif<span class="hljs-punctuation">(</span><span class="hljs-number">500</span><span class="hljs-punctuation">,</span> <span class="hljs-number">0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 标记显著基因</span>
de <span class="hljs-operator">&#x3C;-</span> de <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>significance <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">abs</span><span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span> <span class="hljs-operator">></span> <span class="hljs-number">1</span> <span class="hljs-operator">&#x26;</span> pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"up/down"</span><span class="hljs-punctuation">,</span>
    pvalue <span class="hljs-operator">&#x3C;</span> <span class="hljs-number">0.05</span> <span class="hljs-operator">~</span> <span class="hljs-string">"significant"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span> <span class="hljs-operator">~</span> <span class="hljs-string">"ns"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

ggplot<span class="hljs-punctuation">(</span>de<span class="hljs-punctuation">,</span> aes<span class="hljs-punctuation">(</span>x <span class="hljs-operator">=</span> log2fc<span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span>pvalue<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> significance<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_point<span class="hljs-punctuation">(</span>size <span class="hljs-operator">=</span> <span class="hljs-number">1.5</span><span class="hljs-punctuation">,</span> alpha <span class="hljs-operator">=</span> <span class="hljs-number">0.6</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  scale_color_manual<span class="hljs-punctuation">(</span>values <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"up/down"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#ef4444"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"significant"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#3b82f6"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"ns"</span> <span class="hljs-operator">=</span> <span class="hljs-string">"#9ca3af"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_vline<span class="hljs-punctuation">(</span>xintercept <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-operator">-</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  geom_hline<span class="hljs-punctuation">(</span>yintercept <span class="hljs-operator">=</span> <span class="hljs-operator">-</span>log10<span class="hljs-punctuation">(</span><span class="hljs-number">0.05</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> linetype <span class="hljs-operator">=</span> <span class="hljs-string">"dashed"</span><span class="hljs-punctuation">,</span> color <span class="hljs-operator">=</span> <span class="hljs-string">"gray50"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  labs<span class="hljs-punctuation">(</span>title <span class="hljs-operator">=</span> <span class="hljs-string">"差异表达火山图"</span><span class="hljs-punctuation">,</span> x <span class="hljs-operator">=</span> <span class="hljs-string">"log2 Fold Change"</span><span class="hljs-punctuation">,</span> y <span class="hljs-operator">=</span> <span class="hljs-string">"-log10(p-value)"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">+</span>
  theme_minimal<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>7. 小结</h2>
<ul>
<li>图层语法：<code>ggplot(data, aes()) + geom_xxx() + labs() + theme_xxx()</code></li>
<li>常用几何：<code>geom_point</code> / <code>geom_boxplot</code> / <code>geom_histogram</code> / <code>geom_col</code> / <code>geom_line</code></li>
<li>分组用 <code>aes(color/fill)</code>，拆图用 <code>facet_wrap</code> / <code>facet_grid</code></li>
<li>主题 <code>theme_minimal</code> + <code>labs</code> + <code>ggsave(dpi=300)</code> 即出版级输出</li>
</ul>
<p>至此 R 教程三部曲完成：入门 → tidyverse → ggplot2。</p>`,L_=`<h1>R tidyverse Data Manipulation: dplyr Pipelines and Data Cleaning</h1>
<p>tidyverse is a collection of R packages with a consistent design philosophy, where <code>dplyr</code> handles data manipulation and <code>tidyr</code> handles data reshaping. The core idea is <strong>chaining verbs with pipes</strong>, making every step of data processing clear and readable.</p>
<h2>1. Installation and Loading</h2>
<pre><code class="hljs language-r">install.packages<span class="hljs-punctuation">(</span><span class="hljs-string">"tidyverse"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># install the whole suite at once</span>

library<span class="hljs-punctuation">(</span>tidyverse<span class="hljs-punctuation">)</span>
<span class="hljs-comment"># After loading, automatically includes: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats, tibble</span>
</code></pre>
<h2>2. The Pipe Operator %>%</h2>
<p>The pipe passes the result on the left as the first argument of the function on the right, turning "nested calls" into a "pipeline":</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># Nested style (hard to read)</span>
<span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span>mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Pipe style (read from left to right)</span>
<span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<ul>
<li>Both <code>%>%</code> (magrittr) and <code>|></code> (native, R 4.1+) are available; the tidyverse ecosystem conventionally uses <code>%>%</code></li>
<li>Shortcut: In RStudio, press <code>Ctrl+Shift+M</code> to insert <code>%>%</code></li>
</ul>
<h2>3. Data Preparation</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># Create example data: gene expression experiment</span>
expr_data <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene      <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> each <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  condition <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> length.out <span class="hljs-operator">=</span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  replicate <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">11.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span>    <span class="hljs-comment"># TP53 control</span>
                 <span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">13.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">12.1</span><span class="hljs-punctuation">,</span>   <span class="hljs-comment"># TP53 treatment</span>
                 <span class="hljs-number">8.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># BRCA1</span>
                 <span class="hljs-number">45.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">44.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.8</span><span class="hljs-punctuation">,</span>
                 <span class="hljs-number">33.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">32.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">33.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>expr_data<span class="hljs-punctuation">)</span>
</code></pre>
<p><code>tibble</code> is tidyverse's enhanced data frame: prints more friendlily, and does not automatically convert strings to factors.</p>
<h2>4. Core Verbs</h2>
<h3>4.1 filter: Select Rows</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Samples with expression greater than 20</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">20</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Multiple conditions: &#x26; for AND, | for OR</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span> <span class="hljs-operator">&#x26;</span> <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># %in% for matching a set</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>gene <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 select: Choose Columns</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Exclude columns</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span><span class="hljs-operator">-</span>replicate<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Helper selectors</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"expr"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>      <span class="hljs-comment"># column names starting with "expr"</span>
</code></pre>
<h3>4.3 mutate: Add/Modify Columns</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>
    log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    group <span class="hljs-operator">=</span> paste<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">,</span> sep <span class="hljs-operator">=</span> <span class="hljs-string">"_"</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Conditional column generation: case_when</span>
expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>level <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">30</span> <span class="hljs-operator">~</span> <span class="hljs-string">"high"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span> <span class="hljs-operator">~</span> <span class="hljs-string">"medium"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span>            <span class="hljs-operator">~</span> <span class="hljs-string">"low"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.4 arrange: Sort Rows</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># ascending</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>        <span class="hljs-comment"># descending</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># multi-level sorting</span>
</code></pre>
<h3>4.5 summarise: Summary Statistics</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.6 group_by + summarise: Grouped Summaries (Most Common Combination)</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.7 A Full Pipe Chain</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># Average expression per gene under treatment, top 3 in descending order</span>
top_genes <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  slice_max<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">,</span> n <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>top_genes<span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. Joining Multiple Tables</h2>
<pre><code class="hljs language-r">gene_info <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  function_desc <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"tumor suppressor"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"DNA repair"</span><span class="hljs-punctuation">,</span>
                    <span class="hljs-string">"receptor kinase"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"transcription factor"</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Left join: keeps all rows from the left table</span>
expr_data <span class="hljs-operator">%>%</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Inner join: keeps only rows found in both tables</span>
expr_data <span class="hljs-operator">%>%</span>
  inner_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># Right join / Full join</span>
expr_data <span class="hljs-operator">%>%</span> right_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
expr_data <span class="hljs-operator">%>%</span> full_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>6. tidyr: Data Reshaping</h2>
<h3>6.1 pivot_longer: Wide to Long</h3>
<pre><code class="hljs language-r">wide <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample1 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample2 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample3 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

long <span class="hljs-operator">&#x3C;-</span> wide <span class="hljs-operator">%>%</span>
  pivot_longer<span class="hljs-punctuation">(</span>cols <span class="hljs-operator">=</span> starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"sample"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
               names_to <span class="hljs-operator">=</span> <span class="hljs-string">"sample"</span><span class="hljs-punctuation">,</span>
               values_to <span class="hljs-operator">=</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>long<span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.2 pivot_wider: Long to Wide</h3>
<pre><code class="hljs language-r">long <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> sample<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.3 Handling Missing Values</h3>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> b <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

df <span class="hljs-operator">%>%</span> drop_na<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># remove rows with missing values</span>
df <span class="hljs-operator">%>%</span> replace_na<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-number">0</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># fill specified columns</span>
df <span class="hljs-operator">%>%</span> fill<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>                <span class="hljs-comment"># fill downward</span>
</code></pre>
<h2>7. Practical Example: A Complete Differential Analysis Workflow</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 1. Read data</span>
<span class="hljs-comment"># expr_data &#x3C;- read_csv("expression.csv")</span>

<span class="hljs-comment"># 2. Data cleaning</span>
cleaned <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>              <span class="hljs-comment"># remove missing values</span>
  mutate<span class="hljs-punctuation">(</span>log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>     <span class="hljs-comment"># transformation</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># annotation</span>

<span class="hljs-comment"># 3. Summary statistics</span>
summary_table <span class="hljs-operator">&#x3C;-</span> cleaned <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_log2 <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>log2expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 4. Calculate differences between treatment and control</span>
diff_table <span class="hljs-operator">&#x3C;-</span> summary_table <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> condition<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> mean_log2<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>log2fc <span class="hljs-operator">=</span> treatment <span class="hljs-operator">-</span> control<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>diff_table<span class="hljs-punctuation">)</span>
</code></pre>
<h2>8. Summary</h2>
<ul>
<li><strong>Pipes</strong>: <code>%>%</code> lets data flow from left to right, one verb per step</li>
<li><strong>The six key dplyr verbs</strong>: <code>filter</code> / <code>select</code> / <code>mutate</code> / <code>arrange</code> / <code>summarise</code> / <code>group_by</code></li>
<li><strong>Join family</strong>: <code>left_join</code> / <code>inner_join</code>, etc., merge tables by keys</li>
<li><strong>tidyr</strong>: <code>pivot_longer</code> / <code>pivot_wider</code> convert between wide and long tables</li>
</ul>
<p>The next article will introduce ggplot2: creating publication-quality graphics using the grammar of layers.</p>`,M_=`<h1>R tidyverse 数据操作：dplyr 管道与数据清洗</h1>
<p>tidyverse 是一套风格统一的 R 包集合，其中 <code>dplyr</code> 负责数据操作、<code>tidyr</code> 负责数据整形。核心思想是<strong>管道串联动词</strong>，让数据处理的每一步都清晰可读。</p>
<h2>1. 安装与加载</h2>
<pre><code class="hljs language-r">install.packages<span class="hljs-punctuation">(</span><span class="hljs-string">"tidyverse"</span><span class="hljs-punctuation">)</span>   <span class="hljs-comment"># 一次安装全家桶</span>

library<span class="hljs-punctuation">(</span>tidyverse<span class="hljs-punctuation">)</span>
<span class="hljs-comment"># 加载后自动附带: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats, tibble</span>
</code></pre>
<h2>2. 管道符 %>%</h2>
<p>管道把左边的结果作为右边函数的第一个参数，把"嵌套调用"变成"流水线"：</p>
<pre><code class="hljs language-r"><span class="hljs-comment"># 嵌套写法（难读）</span>
<span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span>mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 管道写法（从左到右读）</span>
<span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">sqrt</span><span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> mean<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">|></span> <span class="hljs-built_in">round</span><span class="hljs-punctuation">(</span><span class="hljs-number">2</span><span class="hljs-punctuation">)</span>
</code></pre>
<ul>
<li><code>%>%</code>（magrittr）与 <code>|></code>（R 4.1+ 原生）都可用；tidyverse 生态惯用 <code>%>%</code></li>
<li>快捷键：RStudio 中 <code>Ctrl+Shift+M</code> 插入 <code>%>%</code></li>
</ul>
<h2>3. 数据准备</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 创建示例数据：基因表达实验</span>
expr_data <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene      <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> each <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  condition <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"control"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> length.out <span class="hljs-operator">=</span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  replicate <span class="hljs-operator">=</span> <span class="hljs-built_in">rep</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-operator">:</span><span class="hljs-number">3</span><span class="hljs-punctuation">,</span> <span class="hljs-number">4</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  <span class="hljs-built_in">expression</span> <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">11.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span>    <span class="hljs-comment"># TP53 control</span>
                 <span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">13.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">12.1</span><span class="hljs-punctuation">,</span>   <span class="hljs-comment"># TP53 treatment</span>
                 <span class="hljs-number">8.1</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">,</span>      <span class="hljs-comment"># BRCA1</span>
                 <span class="hljs-number">45.0</span><span class="hljs-punctuation">,</span> <span class="hljs-number">44.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">45.8</span><span class="hljs-punctuation">,</span>
                 <span class="hljs-number">33.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">32.9</span><span class="hljs-punctuation">,</span> <span class="hljs-number">33.1</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>expr_data<span class="hljs-punctuation">)</span>
</code></pre>
<p><code>tibble</code> 是 tidyverse 的增强数据框：打印更友好、不自动转字符串为因子。</p>
<h2>4. 核心动词</h2>
<h3>4.1 filter：筛选行</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 表达量大于 20 的样本</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">20</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 多条件：&#x26; 且，| 或</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span> <span class="hljs-operator">&#x26;</span> <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">12</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># %in% 匹配集合</span>
expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>gene <span class="hljs-operator">%in%</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.2 select：选择列</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 排除列</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span><span class="hljs-operator">-</span>replicate<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 辅助选择器</span>
expr_data <span class="hljs-operator">%>%</span>
  select<span class="hljs-punctuation">(</span>starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"expr"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>      <span class="hljs-comment"># 列名以 expr 开头</span>
</code></pre>
<h3>4.3 mutate：新增/修改列</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>
    log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    group <span class="hljs-operator">=</span> paste<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">,</span> sep <span class="hljs-operator">=</span> <span class="hljs-string">"_"</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 条件生成列：case_when</span>
expr_data <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>level <span class="hljs-operator">=</span> case_when<span class="hljs-punctuation">(</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">30</span> <span class="hljs-operator">~</span> <span class="hljs-string">"high"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-built_in">expression</span> <span class="hljs-operator">></span> <span class="hljs-number">10</span> <span class="hljs-operator">~</span> <span class="hljs-string">"medium"</span><span class="hljs-punctuation">,</span>
    <span class="hljs-literal">TRUE</span>            <span class="hljs-operator">~</span> <span class="hljs-string">"low"</span>
  <span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.4 arrange：排序</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># 升序</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>        <span class="hljs-comment"># 降序</span>

expr_data <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> desc<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># 多级排序</span>
</code></pre>
<h3>4.5 summarise：汇总</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.6 group_by + summarise：分组汇总（最常用组合）</h3>
<pre><code class="hljs language-r">expr_data <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    sd_expr <span class="hljs-operator">=</span> sd<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>
</code></pre>
<h3>4.7 完整管道链</h3>
<pre><code class="hljs language-r"><span class="hljs-comment"># 每个基因在 treatment 下的平均表达，按降序取前 3</span>
top_genes <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span>condition <span class="hljs-operator">==</span> <span class="hljs-string">"treatment"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>mean_expr <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  slice_max<span class="hljs-punctuation">(</span>mean_expr<span class="hljs-punctuation">,</span> n <span class="hljs-operator">=</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>top_genes<span class="hljs-punctuation">)</span>
</code></pre>
<h2>5. 连接多个表（join）</h2>
<pre><code class="hljs language-r">gene_info <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"EGFR"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"MYC"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  chromosome <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">17</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  function_desc <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"tumor suppressor"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"DNA repair"</span><span class="hljs-punctuation">,</span>
                    <span class="hljs-string">"receptor kinase"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"transcription factor"</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 左连接：保留左表所有行</span>
expr_data <span class="hljs-operator">%>%</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 内连接：只保留两边都有的</span>
expr_data <span class="hljs-operator">%>%</span>
  inner_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 右连接 / 全连接</span>
expr_data <span class="hljs-operator">%>%</span> right_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
expr_data <span class="hljs-operator">%>%</span> full_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>
</code></pre>
<h2>6. tidyr：数据整形</h2>
<h3>6.1 pivot_longer：宽表转长表</h3>
<pre><code class="hljs language-r">wide <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>
  gene <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-string">"TP53"</span><span class="hljs-punctuation">,</span> <span class="hljs-string">"BRCA1"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample1 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">10.2</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.1</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample2 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">12.5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">7.9</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
  sample3 <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">9.8</span><span class="hljs-punctuation">,</span> <span class="hljs-number">8.4</span><span class="hljs-punctuation">)</span>
<span class="hljs-punctuation">)</span>

long <span class="hljs-operator">&#x3C;-</span> wide <span class="hljs-operator">%>%</span>
  pivot_longer<span class="hljs-punctuation">(</span>cols <span class="hljs-operator">=</span> starts_with<span class="hljs-punctuation">(</span><span class="hljs-string">"sample"</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
               names_to <span class="hljs-operator">=</span> <span class="hljs-string">"sample"</span><span class="hljs-punctuation">,</span>
               values_to <span class="hljs-operator">=</span> <span class="hljs-string">"expression"</span><span class="hljs-punctuation">)</span>
print<span class="hljs-punctuation">(</span>long<span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.2 pivot_wider：长表转宽表</h3>
<pre><code class="hljs language-r">long <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> sample<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> <span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span>
</code></pre>
<h3>6.3 缺失值处理</h3>
<pre><code class="hljs language-r">df <span class="hljs-operator">&#x3C;-</span> tibble<span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-number">1</span><span class="hljs-punctuation">,</span> <span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">3</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span> b <span class="hljs-operator">=</span> <span class="hljs-built_in">c</span><span class="hljs-punctuation">(</span><span class="hljs-literal">NA</span><span class="hljs-punctuation">,</span> <span class="hljs-number">5</span><span class="hljs-punctuation">,</span> <span class="hljs-number">6</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

df <span class="hljs-operator">%>%</span> drop_na<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span>              <span class="hljs-comment"># 删除含缺失的行</span>
df <span class="hljs-operator">%>%</span> replace_na<span class="hljs-punctuation">(</span><span class="hljs-built_in">list</span><span class="hljs-punctuation">(</span>a <span class="hljs-operator">=</span> <span class="hljs-number">0</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>  <span class="hljs-comment"># 指定列填充</span>
df <span class="hljs-operator">%>%</span> fill<span class="hljs-punctuation">(</span>a<span class="hljs-punctuation">)</span>                <span class="hljs-comment"># 向下填充</span>
</code></pre>
<h2>7. 实战：完整差异分析工作流</h2>
<pre><code class="hljs language-r"><span class="hljs-comment"># 1. 读取数据</span>
<span class="hljs-comment"># expr_data &#x3C;- read_csv("expression.csv")</span>

<span class="hljs-comment"># 2. 数据清洗</span>
cleaned <span class="hljs-operator">&#x3C;-</span> expr_data <span class="hljs-operator">%>%</span>
  filter<span class="hljs-punctuation">(</span><span class="hljs-operator">!</span><span class="hljs-built_in">is.na</span><span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>              <span class="hljs-comment"># 去缺失</span>
  mutate<span class="hljs-punctuation">(</span>log2expr <span class="hljs-operator">=</span> log2<span class="hljs-punctuation">(</span><span class="hljs-built_in">expression</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>     <span class="hljs-comment"># 转换</span>
  left_join<span class="hljs-punctuation">(</span>gene_info<span class="hljs-punctuation">,</span> by <span class="hljs-operator">=</span> <span class="hljs-string">"gene"</span><span class="hljs-punctuation">)</span>           <span class="hljs-comment"># 注释</span>

<span class="hljs-comment"># 3. 汇总统计</span>
summary_table <span class="hljs-operator">&#x3C;-</span> cleaned <span class="hljs-operator">%>%</span>
  group_by<span class="hljs-punctuation">(</span>gene<span class="hljs-punctuation">,</span> condition<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  summarise<span class="hljs-punctuation">(</span>
    mean_log2 <span class="hljs-operator">=</span> mean<span class="hljs-punctuation">(</span>log2expr<span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    n <span class="hljs-operator">=</span> n<span class="hljs-punctuation">(</span><span class="hljs-punctuation">)</span><span class="hljs-punctuation">,</span>
    .groups <span class="hljs-operator">=</span> <span class="hljs-string">"drop"</span>
  <span class="hljs-punctuation">)</span>

<span class="hljs-comment"># 4. 计算处理/对照的差异</span>
diff_table <span class="hljs-operator">&#x3C;-</span> summary_table <span class="hljs-operator">%>%</span>
  pivot_wider<span class="hljs-punctuation">(</span>names_from <span class="hljs-operator">=</span> condition<span class="hljs-punctuation">,</span> values_from <span class="hljs-operator">=</span> mean_log2<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  mutate<span class="hljs-punctuation">(</span>log2fc <span class="hljs-operator">=</span> treatment <span class="hljs-operator">-</span> control<span class="hljs-punctuation">)</span> <span class="hljs-operator">%>%</span>
  arrange<span class="hljs-punctuation">(</span>desc<span class="hljs-punctuation">(</span>log2fc<span class="hljs-punctuation">)</span><span class="hljs-punctuation">)</span>

print<span class="hljs-punctuation">(</span>diff_table<span class="hljs-punctuation">)</span>
</code></pre>
<h2>8. 小结</h2>
<ul>
<li><strong>管道</strong>：<code>%>%</code> 让数据流从左到右，每一步一个动词</li>
<li><strong>dplyr 六大动词</strong>：<code>filter</code> / <code>select</code> / <code>mutate</code> / <code>arrange</code> / <code>summarise</code> / <code>group_by</code></li>
<li><strong>join 家族</strong>：<code>left_join</code> / <code>inner_join</code> 等按键合并表</li>
<li><strong>tidyr</strong>：<code>pivot_longer</code> / <code>pivot_wider</code> 在宽长表之间转换</li>
</ul>
<p>下一篇将介绍 ggplot2：用图层语法绘制出版级图表。</p>`,D_=`<h1>Review of Cryo-EM Technology</h1>
<p>In 2015, Nature Methods named cryo-electron microscopy (cryo-EM) its "Method of the Year." A decade on, this technique has evolved from the subject of ridicule as "blobology" into a mainstay for determining the structures of biological macromolecules. This article reviews the key milestones of this technological revolution and looks ahead to structural biology in the AI era.</p>
<h2>1. The Eve of Revolution: Why Cryo-EM Was Once Called "Blobology"</h2>
<h3>1.1 The Three Pillars of Traditional Structural Biology</h3>
<ul>
<li><strong>X-ray crystallography</strong>: Requires high-quality crystals; many membrane proteins and large complexes are difficult to crystallize</li>
<li><strong>Nuclear magnetic resonance (NMR)</strong>: Limited by molecular weight (typically &#x3C; 50 kDa)</li>
<li><strong>Cryo-EM</strong>: No upper molecular weight limit in theory, no crystals required, but resolution remained stuck at 10–30 Å for a long time</li>
</ul>
<h3>1.2 Root Causes of the Resolution Bottleneck</h3>
<ol>
<li><strong>Electron damage</strong>: The electron beam damages the sample, limiting the total dose to ~20 e⁻/Å², resulting in extremely low signal-to-noise ratio</li>
<li><strong>Outdated detectors</strong>: CCD cameras are noisy at low doses and lack single-electron counting capability</li>
<li><strong>Specimen drift</strong>: Electron-beam-induced specimen movement blurs the images</li>
<li><strong>Missing phase information</strong>: Under the weak-phase-object approximation, defocused imaging loses part of the phase information</li>
</ol>
<h2>2. The Three Cornerstones of the Technological Revolution</h2>
<h3>2.1 Direct Electron Detectors (DET) — 2013–2015</h3>
<table>
<thead>
<tr>
<th>Detector</th>
<th>Year Released</th>
<th>Key Features</th>
</tr>
</thead>
<tbody>
<tr>
<td>Gatan K2 Summit</td>
<td>2013</td>
<td>First commercial CMOS direct detector, single-electron counting</td>
</tr>
<tr>
<td>FEI Falcon II</td>
<td>2014</td>
<td>High frame rate, lossless readout</td>
</tr>
<tr>
<td>Gatan K3</td>
<td>2017</td>
<td>Larger field of view, super-resolution</td>
</tr>
<tr>
<td>Falcon 4 / Falcon 4i</td>
<td>2020</td>
<td>Electron counting + integrated energy filtering</td>
</tr>
</tbody>
</table>
<p><strong>Why DET is revolutionary</strong>:</p>
<ul>
<li>Single-electron counting mode eliminates readout noise → greatly improved signal-to-noise ratio</li>
<li>Ultra-fast frame rates (hundreds of frames per second) → <strong>movie mode</strong>, enabling per-frame correction of specimen drift</li>
<li>Combined with energy filters (GIF/Selectris), removes inelastic scattering background</li>
</ul>
<h3>2.2 The Stability Revolution — 2016–2018</h3>
<p>The "Resolution Revolution" was not just about detectors; it was equally a triumph of <strong>overall stability</strong>:</p>
<ul>
<li>Constant temperature and humidity in the microscope room (±0.1 °C), soundproofing and vibration isolation</li>
<li>Improved specimen holders and grid-clamping mechanisms (e.g., autoloader)</li>
<li>Widespread adoption of dose-symmetric acquisition schemes and low-dose modes</li>
<li>Software: MotionCor2's <strong>per-frame drift correction</strong> preserved high-frequency signals</li>
</ul>
<p>Representative achievements (2017–2020):</p>
<ul>
<li>Human γ-secretase complex (3.4 Å, 2015, Scheres group)</li>
<li>β-galactosidase at 2.2 Å (2016)</li>
<li><strong>Human ribosome at ~2.5 Å</strong> (multiple groups)</li>
<li>Multiple membrane protein channels and receptor complexes broke the 3 Å barrier</li>
</ul>
<h3>2.3 Algorithms and Computing Power — 2017–2021</h3>
<ul>
<li><strong>Bayesian methods</strong> (RELION) made 3D classification/refinement robust</li>
<li><strong>cryoSPARC's ab-initio and NU-refinement</strong>: heterogeneous species separation, flexible-region refinement</li>
<li><strong>GPU acceleration</strong> made iterative reconstruction of million-particle datasets routine</li>
<li>Deep learning enters: <strong>Topaz / crYOLO</strong> for particle picking, <strong>DeepEM</strong> for denoising</li>
</ul>
<h2>3. The 2020s: From 3 Å to Atomic Resolution</h2>
<h3>3.1 Breaking the 2 Å Barrier</h3>
<ul>
<li>2020–2021: Multiple samples reached 1.5–2.0 Å (e.g., apoferritin at 1.22 Å, 2020)</li>
<li>Sub-2 Å density maps reveal: hydrogen atoms, water molecule networks, side-chain microenvironments, ligand interaction details</li>
<li>This paves the way for <strong>density-based de novo modeling</strong> and drug design</li>
</ul>
<h3>3.2 Time-Resolved Cryo-EM</h3>
<ul>
<li>Mix-spray techniques capture enzyme catalytic intermediates</li>
<li>Microfluidics and rapid freezing (&#x3C; 10 ms timescale)</li>
<li>Goal: observe conformational change dynamics, rather than a single static structure</li>
</ul>
<h3>3.3 In Situ Structural Biology (In situ / Cellular Cryo-EM)</h3>
<ul>
<li>Cryo-focused ion beam (cryo-FIB) thinning + cryo-electron tomography (cryo-ET)</li>
<li>Resolving the assembly states and spatial distribution of protein complexes <strong>within cells</strong></li>
<li>Representative: in situ structures of ribosomes, proteasomes, and vesicle trafficking complexes in cells</li>
</ul>
<h2>4. Synergy and Challenges in the AI Era</h2>
<h3>4.1 The AlphaFold Shockwave (2021–)</h3>
<p>AlphaFold2 / AlphaFold3 can predict protein structures with high accuracy, forming a <strong>complementary rather than replacement</strong> relationship with cryo-EM:</p>
<table>
<thead>
<tr>
<th>Dimension</th>
<th>AlphaFold</th>
<th>Cryo-EM</th>
</tr>
</thead>
<tbody>
<tr>
<td>Input</td>
<td>Sequence</td>
<td>Experimental sample</td>
</tr>
<tr>
<td>Output</td>
<td>Predicted model</td>
<td>Experimental density/model</td>
</tr>
<tr>
<td>Conformation</td>
<td>Single (or few)</td>
<td>True conformational distribution</td>
</tr>
<tr>
<td>Modifications/ligands</td>
<td>Limited predictive capability</td>
<td>Direct observation</td>
</tr>
<tr>
<td>Complex assembly</td>
<td>Predicted</td>
<td>True in situ state</td>
</tr>
</tbody>
</table>
<p><strong>Actual workflow</strong>: AF2 predicted model → used as initial template for molecular replacement/model building in cryo-EM → refined against experimental density → yielding a model containing true interactions and modifications.</p>
<h3>4.2 Further Penetration of Deep Learning</h3>
<ul>
<li>Automatic density map model building (ModelAngelo, DeepTracer)</li>
<li>End-to-end models for particle picking and denoising</li>
<li>Diffusion-model-based density map deblurring</li>
<li><strong>Joint refinement of AlphaFold and density maps</strong> (e.g., ISOLDE + AF constraints)</li>
</ul>
<h3>4.3 Challenges Remain</h3>
<ul>
<li><strong>Small proteins (&#x3C; 50 kDa)</strong>: low signal-to-noise ratio, difficult to improve resolution (require better phase plates/detectors)</li>
<li><strong>Flexible large complexes</strong>: conformational continua are difficult to describe with discrete classification (3DVA, 3D FLEX, etc. are addressing this)</li>
<li><strong>Membrane proteins</strong>: detergent/nanodisc environments differ from native membranes</li>
<li><strong>Specimen preparation</strong>: the air-water interface problem remains one of the major bottlenecks</li>
</ul>
<h2>5. Future Outlook</h2>
<ol>
<li><strong>Resolution becoming routine</strong>: sub-2 Å transitions from "world record" to a routine goal</li>
<li><strong>Dynamic structural biology</strong>: time-resolved + continuous conformational analysis, moving from "snapshots" to "movies"</li>
<li><strong>In-cell structural biology</strong>: cryo-ET combined with FIB to understand molecular machines in the cellular environment</li>
<li><strong>Full AI integration</strong>: automated, intelligent pipelines from picking to modeling</li>
<li><strong>Deep integration with drug design</strong>: high-resolution structures + virtual screening + generative design closed loop</li>
</ol>
<h2>6. Summary</h2>
<ul>
<li>Direct electron detectors, the stability revolution, and algorithmic advances are the three cornerstones of the "Resolution Revolution"</li>
<li>Cryo-EM has transformed from a "last resort" into the <strong>preferred method</strong> for membrane proteins and large complexes</li>
<li>AlphaFold and cryo-EM are complementary: prediction guides experiment, experiment corrects prediction</li>
<li>Future directions: dynamics, in situ, atomic resolution, AI integration</li>
</ul>
<p>The past decade of cryo-EM is a model of synergy among engineering, physics, and computational biology. In the next decade, it will join forces with AI to redefine how we understand life.</p>`,O_=`<h1>冷冻电镜技术综述</h1>
<p>2015 年，Nature Methods 将冷冻电镜（cryo-EM）评为"年度方法"。十年过去，这项技术从"blobology"（模糊团块学）的调侃对象，成长为解析生物大分子结构的主力军。本文回顾这场技术革命的关键节点，并展望 AI 时代的结构生物学。</p>
<h2>1. 革命前夜：为什么冷冻电镜曾被称为"blobology"</h2>
<h3>1.1 传统结构生物学的三驾马车</h3>
<ul>
<li><strong>X 射线晶体学</strong>：需要高质量晶体，许多膜蛋白、大复合物难以结晶</li>
<li><strong>核磁共振（NMR）</strong>：受分子量限制（通常 &#x3C; 50 kDa）</li>
<li><strong>冷冻电镜</strong>：理论上无分子量上限、无需晶体，但分辨率长期卡在 10–30 Å</li>
</ul>
<h3>1.2 分辨率瓶颈的根源</h3>
<ol>
<li><strong>电子损伤</strong>：电子束会破坏样品，总剂量被限制在 ~20 e⁻/Å²，信噪比极低</li>
<li><strong>探测器落后</strong>：CCD 相机在低剂量下噪声大，且无单电子计数能力</li>
<li><strong>样品漂移</strong>：电子束导致的样品移动模糊图像</li>
<li><strong>相位信息缺失</strong>：弱相位物体近似下，欠焦成像丢失部分相位信息</li>
</ol>
<h2>2. 技术革命的三块基石</h2>
<h3>2.1 直接电子探测器（DET）——2013–2015</h3>
<table>
<thead>
<tr>
<th>探测器</th>
<th>推出年份</th>
<th>关键特性</th>
</tr>
</thead>
<tbody>
<tr>
<td>Gatan K2 Summit</td>
<td>2013</td>
<td>首款商用 CMOS 直接探测，单电子计数</td>
</tr>
<tr>
<td>FEI Falcon II</td>
<td>2014</td>
<td>高帧率、无损读出</td>
</tr>
<tr>
<td>Gatan K3</td>
<td>2017</td>
<td>更大视野、超分辨率</td>
</tr>
<tr>
<td>Falcon 4 / Falcon 4i</td>
<td>2020</td>
<td>电子计数 + 能量过滤集成</td>
</tr>
</tbody>
</table>
<p><strong>为什么 DET 是革命性的</strong>：</p>
<ul>
<li>单电子计数模式消除读出噪声 → 信噪比大幅提升</li>
<li>超快帧率（每秒数百帧）→ <strong>movie 模式</strong>，可逐帧校正样品漂移</li>
<li>与能量过滤器（GIF/Selectris）结合，去除非弹性散射本底</li>
</ul>
<h3>2.2 稳定性革命——2016–2018</h3>
<p>"Resolution Revolution" 不只是探测器，更是<strong>整体稳定性</strong>的胜利：</p>
<ul>
<li>电镜室恒温恒湿（±0.1 °C）、隔音防震</li>
<li>样品杆与载网夹持机构改进（如 autoloader）</li>
<li>剂量对称采集方案、低剂量模式的普及</li>
<li>软件层面：MotionCor2 的<strong>逐帧漂移校正</strong>让高频信号得以保留</li>
</ul>
<p>代表性成果（2017–2020）：</p>
<ul>
<li>人源 γ-分泌酶复合物（3.4 Å，2015，Scheres 组）</li>
<li>β-半乳糖苷酶 2.2 Å（2016）</li>
<li><strong>人源核糖体 2.5 Å 级</strong>（多组）</li>
<li>多种膜蛋白通道、受体复合物突破 3 Å</li>
</ul>
<h3>2.3 算法与算力——2017–2021</h3>
<ul>
<li><strong>贝叶斯方法</strong>（RELION）让 3D 分类/精修稳健化</li>
<li><strong>cryoSPARC 的 Ab-initio 与 NU-refinement</strong>：异构体分离、柔性区域精修</li>
<li><strong>GPU 加速</strong>使百万级颗粒的迭代重构成为常规操作</li>
<li>深度学习进入：<strong>Topaz / crYOLO</strong> 颗粒挑选、<strong>DeepEM</strong> 去噪</li>
</ul>
<h2>3. 2020 年代：从 3 Å 到原子级</h2>
<h3>3.1 打破 2 Å 门槛</h3>
<ul>
<li>2020–2021：多种样品达到 1.5–2.0 Å（如去铁蛋白 1.22 Å，2020）</li>
<li>亚 2 Å 密度下可观察到：氢原子、水分子网络、侧链微环境、配体互作细节</li>
<li>这为<strong>基于密度的从头建模</strong>与药物设计铺平道路</li>
</ul>
<h3>3.2 时间分辨冷冻电镜（Time-resolved cryo-EM）</h3>
<ul>
<li>混合-喷雾（mix-spray）技术捕捉酶催化中间态</li>
<li>微流控与快速冷冻（&#x3C; 10 ms 时间尺度）</li>
<li>目标：观察构象变化动态过程，而非单一静态结构</li>
</ul>
<h3>3.3 原位结构生物学（In situ / Cellular cryo-EM）</h3>
<ul>
<li>冷冻聚焦离子束（cryo-FIB）减薄 + 冷冻电子断层扫描（cryo-ET）</li>
<li>在<strong>细胞原位</strong>解析蛋白质复合物的组装状态与空间分布</li>
<li>代表：细胞内核糖体、蛋白酶体、囊泡运输复合物的原位结构</li>
</ul>
<h2>4. AI 时代的协同与挑战</h2>
<h3>4.1 AlphaFold 冲击波（2021–）</h3>
<p>AlphaFold2 / AlphaFold3 可以高精度预测蛋白结构，与冷冻电镜形成<strong>互补而非替代</strong>：</p>
<table>
<thead>
<tr>
<th>维度</th>
<th>AlphaFold</th>
<th>冷冻电镜</th>
</tr>
</thead>
<tbody>
<tr>
<td>输入</td>
<td>序列</td>
<td>实验样品</td>
</tr>
<tr>
<td>输出</td>
<td>预测模型</td>
<td>实验密度/模型</td>
</tr>
<tr>
<td>构象</td>
<td>单一（或少量）</td>
<td>真实构象分布</td>
</tr>
<tr>
<td>修饰/配体</td>
<td>预测能力有限</td>
<td>直接观察</td>
</tr>
<tr>
<td>复合物组装</td>
<td>预测</td>
<td>原位真实状态</td>
</tr>
</tbody>
</table>
<p><strong>实际工作流</strong>：AF2 预测模型 → 作为冷冻电镜分子替换（molecular replacement）/建模的初始模板 → 实验密度修正 → 获得包含真实互作与修饰的模型。</p>
<h3>4.2 深度学习的进一步渗透</h3>
<ul>
<li>密度图自动建模（ModelAngelo、DeepTracer）</li>
<li>颗粒挑选与去噪的端到端模型</li>
<li>基于扩散模型的密度图去模糊</li>
<li><strong>AlphaFold 与密度图的联合精修</strong>（如 ISOLDE + AF 约束）</li>
</ul>
<h3>4.3 挑战依然存在</h3>
<ul>
<li><strong>小蛋白（&#x3C; 50 kDa）</strong>：信噪比低，分辨率提升困难（需要更好的相板/探测器）</li>
<li><strong>柔性大复合物</strong>：构象连续体难以用离散分类描述（3DVA、3D FLEX 等正在解决）</li>
<li><strong>膜蛋白</strong>：去垢剂/纳米盘环境与天然膜差异</li>
<li><strong>样品制备</strong>：air-water interface 问题仍是主要瓶颈之一</li>
</ul>
<h2>5. 未来图景</h2>
<ol>
<li><strong>分辨率常态化</strong>：亚 2 Å 从"世界纪录"变为常规目标</li>
<li><strong>动态结构生物学</strong>：时间分辨 + 连续构象分析，从"快照"走向"电影"</li>
<li><strong>细胞内结构生物学</strong>：cryo-ET 与 FIB 结合，在细胞环境中理解分子机器</li>
<li><strong>AI 全流程集成</strong>：从挑选到建模的自动化、智能化管线</li>
<li><strong>与药物设计深度融合</strong>：高分辨结构 + 虚拟筛选 + 生成式设计闭环</li>
</ol>
<h2>6. 小结</h2>
<ul>
<li>直接电子探测器、稳定性革命、算法进步是"分辨率革命"的三块基石</li>
<li>冷冻电镜已从"最后手段"变为膜蛋白、大复合物的<strong>首选方法</strong></li>
<li>AlphaFold 与冷冻电镜形成互补：预测引导实验，实验修正预测</li>
<li>未来方向：动态、原位、原子级、AI 融合</li>
</ul>
<p>冷冻电镜的十年，是工程、物理与计算生物学协同的典范。下一个十年，它将与 AI 一起重新定义我们理解生命的方式。</p>`,I_=`<h1>Cryo-EM Single Particle Analysis: Full Data Processing Workflow</h1>
<p>Cryo-EM single particle analysis (SPA) is a core technique for determining near-atomic resolution structures of biological macromolecules. This article outlines the complete pipeline from micrographs to atomic models, introducing the purpose, common software, and key parameters for each step.</p>
<h2>1. Workflow Overview</h2>
<pre><code>Sample preparation (cryo-grid preparation)
   ↓
Data acquisition (Titan Krios / Glacios, etc.)
   ↓
① Preprocessing: motion correction + CTF estimation
   ↓
② Particle picking
   ↓
③ 2D classification (removing bad particles)
   ↓
④ Initial model (Ab initio)
   ↓
⑤ 3D classification (heterogeneity analysis)
   ↓
⑥ 3D refinement
   ↓
⑦ Resolution assessment (FSC curve)
   ↓
⑧ Atomic model building and refinement
</code></pre>
<h2>2. Sample Preparation and Data Acquisition</h2>
<h3>2.1 Cryo-grid Preparation</h3>
<ul>
<li>Apply the protein sample (typically 0.5–5 mg/mL) to a holey carbon grid</li>
<li>Use a Vitrobot / Leica GP or similar device to rapidly plunge into liquid ethane, forming <strong>vitreous ice</strong></li>
<li>Ideal ice thickness: slightly thinner than the particle diameter; ice that is too thick reduces signal-to-noise, while ice that is too thin may deform particles</li>
</ul>
<h3>2.2 Data Acquisition</h3>
<ul>
<li>300 kV electron microscope (Titan Krios) + direct electron detector (Gatan K3 / Falcon 4)</li>
<li><strong>Super-resolution mode</strong> + dose fractionation (total dose about 40–60 e⁻/Å²)</li>
<li>Each movie typically contains 30–50 frames for motion correction</li>
<li>Target: 10,000–20,000 micrographs (high-resolution structures usually require millions of particles)</li>
</ul>
<h2>3. Preprocessing</h2>
<h3>3.1 Motion Correction</h3>
<p>Sample drift caused by the electron beam blurs images and must be corrected frame by frame:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Motion correction in RELION 4</span>
relion_run_motioncorr \\
  --i Micrographs/ \\
  --o MotionCorr/ \\
  --j 8 \\
  --pipeline_control MotionCorr/job001/
</code></pre>
<p>Common software: <strong>MotionCor2</strong> (classic), <strong>built-in RELION</strong>, <strong>cryoSPARC Patch Motion Correction</strong>.</p>
<h3>3.2 CTF Estimation</h3>
<p>The CTF (Contrast Transfer Function) describes the phase flipping and defocus effects of electron microscope imaging:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># CTFFIND4 (classic command line)</span>
ctffind --input_micrograph MotionCorr/mic1.mrc \\
        --output_diag mic1_diag.mrc \\
        --defU 20000 --defV 20000

<span class="hljs-comment"># Gctf (faster GPU version)</span>
Gctf --input <span class="hljs-string">'MotionCorr/*.mrc'</span> --apix 1.0
</code></pre>
<p>Common software: <strong>CTFFIND4</strong>, <strong>Gctf</strong>, <strong>cryoSPARC Patch CTF</strong>.</p>
<p><strong>Key criteria</strong>: The resolution cutoff is usually taken as 3–4 Å; discard micrographs with excessive astigmatism or severe drift.</p>
<h2>4. Particle Picking</h2>
<h3>4.1 Traditional Methods</h3>
<ul>
<li><strong>Template-based</strong>: use 2D class averages as templates and perform correlation searches</li>
<li>Software: auto-picking in RELION, Template Picker in cryoSPARC</li>
</ul>
<h3>4.2 Deep Learning Approaches (currently mainstream)</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Topaz: CNN-based particle picking</span>
topaz train --model-path model.pkl --epochs 20 \\
  --num-workers 4 --train-data particles.toml

topaz extract --model-path model.pkl \\
  --micrographs micrographs.toml \\
  --output-prefix picked --radius 100
</code></pre>
<p>Mainstream tools: <strong>Topaz</strong>, <strong>crYOLO</strong> (very fast, supports real-time picking), <strong>cryoSPARC Blob Picker</strong>.</p>
<blockquote>
<p>In modern workflows, a blob/template picker is often used first to quickly select a subset for 2D classification, and clean class averages are then used to train Topaz/crYOLO models, significantly improving picking quality for complex samples (membrane proteins, small proteins).</p>
</blockquote>
<h2>5. 2D Classification: Particle Screening</h2>
<p>2D classification clusters particles by projection orientation while <strong>removing ice contamination, mis-picked particles, and denatured particles</strong>:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 2D classification</span>
relion_refine --o Class2D/ \\
  --particle_images particles.star \\
  --ctf --K 100 --iter 25 --tau2_fudge 4 \\
  --flatten_solvent --zero_mask --pool 30 \\
  --per_particle_ctf --j 8
</code></pre>
<p><strong>Quality criteria</strong>:</p>
<ul>
<li>Class averages should show clear secondary structure (α-helices, β-sheets)</li>
<li>Retain 50–80% of the total particles as "good"</li>
<li>High-quality class averages are a prerequisite for subsequent 3D reconstruction</li>
</ul>
<h2>6. 3D Reconstruction: Initial Model and Classification</h2>
<h3>6.1 Initial Model (Ab initio)</h3>
<p>No reference model is required; reconstruction starts from random initial structures:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># cryoSPARC</span>
cryosparc_abinitio ... --num_classes 3
</code></pre>
<ul>
<li>Commonly used: <strong>cryoSPARC Ab-initio</strong> (fast and robust), RELION 3D initial model</li>
<li>Generate 2–4 classes and select the class with the clearest structural features as a reference</li>
</ul>
<h3>6.2 3D Classification: Handling Conformational Heterogeneity</h3>
<p>Biological samples often exist in multiple conformations (open/closed, ligand-bound/unbound):</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 3D classification</span>
relion_refine --o Class3D/ \\
  --ref initial.mrc --K 4 --iter 40 \\
  --tau2_fudge 8 --particle_diameter 200
</code></pre>
<p><strong>Strategy</strong>: First use coarse classification (K=4–6) to inspect the conformational distribution, then sub-classify the conformation of interest (K=2–3). Ultimately, each 3D class corresponds to one conformation.</p>
<h2>7. 3D Refinement</h2>
<h3>7.1 Standard Refinement</h3>
<pre><code class="hljs language-bash">relion_refine --o Refine3D/ \\
  --ref class3d.mrc --particle_diameter 200 \\
  --auto_mask --solvent_mask \\
  --ctf --per_particle_ctf --j 8
</code></pre>
<p>Key parameters:</p>
<ul>
<li>Start with <strong>C1 symmetry</strong>; for symmetric molecules, specify C/D/O symmetry (e.g., C3, D7)</li>
<li>A <strong>solvent mask</strong> must be used; otherwise the FSC is artificially high</li>
<li>Use local refinement to handle flexible regions</li>
</ul>
<h3>7.2 Bayesian Polishing</h3>
<p>Performs trajectory fitting and weighting at the movie-frame level for particles, improving high-frequency information:</p>
<pre><code class="hljs language-bash">relion_polish --i Refine3D/run_data.star \\
  --o Polish/ \\
  --model_type 2 --nr_iter 5
</code></pre>
<h3>7.3 CTF Refinement</h3>
<p>Refines per-particle defocus, astigmatism, and beam tilt:</p>
<pre><code class="hljs language-bash">relion_ctf_refine --i Polish/particles.star \\
  --o CtfRefine/ \\
  --fit_defocus --fit_magnification --fit_tilt
</code></pre>
<h2>8. Resolution Assessment: FSC Curve</h2>
<p><strong>FSC (Fourier Shell Correlation)</strong>: Particles are randomly split into two half-sets (half-maps), reconstructed separately, and the correlation coefficient between the two maps is computed:</p>
<ul>
<li><strong>FSC = 0.143 criterion</strong>: gold-standard resolution determination (after half-map correction)</li>
<li><strong>FSC = 0.5</strong>: conservative estimate (uncorrected)</li>
</ul>
<pre><code class="hljs language-bash">relion_postprocess --i Refine3D/run_half1_class001.mrc \\
  --mask mask.mrc --autob_masking \\
  --angpix 1.0
</code></pre>
<p><strong>Quality checklist</strong>:</p>
<table>
<thead>
<tr>
<th>Metric</th>
<th>Good standard</th>
</tr>
</thead>
<tbody>
<tr>
<td>Global resolution (FSC 0.143)</td>
<td>≤ 3.5 Å (high resolution)</td>
</tr>
<tr>
<td>Number of particles</td>
<td>≥ 100,000 (hundreds of thousands typically needed for 3 Å)</td>
</tr>
<tr>
<td>Orientation distribution</td>
<td>Uniform (no preferred orientation)</td>
</tr>
<tr>
<td>Density continuity</td>
<td>Side-chain density clearly visible (&#x3C; 3 Å)</td>
</tr>
</tbody>
</table>
<h2>9. Comparison of Common Software Stacks</h2>
<table>
<thead>
<tr>
<th>Step</th>
<th>RELION (command line)</th>
<th>cryoSPARC (GUI)</th>
</tr>
</thead>
<tbody>
<tr>
<td>Motion correction</td>
<td>Built-in / MotionCor2</td>
<td>Patch Motion</td>
</tr>
<tr>
<td>CTF</td>
<td>CTFFIND4 / Gctf</td>
<td>Patch CTF</td>
</tr>
<tr>
<td>Picking</td>
<td>Template / Topaz integration</td>
<td>Blob / Template / Topaz</td>
</tr>
<tr>
<td>2D/3D</td>
<td>Stable and reliable</td>
<td>Fast, interactive</td>
</tr>
<tr>
<td>Non-uniform refinement</td>
<td>—</td>
<td><strong>NU-refinement</strong> (powerful for flexible macromolecules)</td>
</tr>
</tbody>
</table>
<blockquote>
<p><strong>Workflow recommendation</strong>: Use cryoSPARC for preprocessing and initial reconstruction, and RELION for high-precision refinement and polishing; alternatively, use cryoSPARC for the entire workflow (more user-friendly for beginners).</p>
</blockquote>
<h2>10. Summary</h2>
<ul>
<li>Preprocessing (motion correction + CTF) determines the upper limit of data quality</li>
<li>Deep learning-based particle picking (Topaz/crYOLO) significantly improves performance on complex samples</li>
<li>2D classification screens particles; 3D classification handles conformational heterogeneity</li>
<li>Three key refinement steps: solvent mask + Bayesian polishing + CTF refinement</li>
<li>FSC = 0.143 is the gold standard for resolution</li>
</ul>
<p>Future articles will cover structure visualization (PyMOL/ChimeraX) and atomic modeling (Coot/phenix), turning density maps into atomic models.</p>`,N_=`<h1>冷冻电镜单颗粒分析：数据处理全流程</h1>
<p>冷冻电镜单颗粒分析（Single Particle Analysis, SPA）是解析生物大分子近原子分辨率结构的核心技术。本文梳理从显微照片到原子模型的完整流程，介绍各步骤的目的、常用软件与关键参数。</p>
<h2>1. 流程总览</h2>
<pre><code>样品制备（冷冻制样）
   ↓
数据采集（Titan Krios / Glacios 等）
   ↓
① 预处理：运动校正 + CTF 估计
   ↓
② 颗粒挑选（Particle Picking）
   ↓
③ 2D 分类（筛除坏颗粒）
   ↓
④ 初始模型（Ab-initio）
   ↓
⑤ 3D 分类（异质性分析）
   ↓
⑥ 3D 精修（Refinement）
   ↓
⑦ 分辨率评估（FSC 曲线）
   ↓
⑧ 原子模型搭建与精修
</code></pre>
<h2>2. 样品制备与数据采集</h2>
<h3>2.1 冷冻制样</h3>
<ul>
<li>将蛋白样品（浓度通常 0.5–5 mg/mL）施加到载网（holey carbon grid）</li>
<li>用 Vitrobot / Leica GP 等设备快速浸入液乙烷，形成<strong>玻璃态冰</strong></li>
<li>理想冰厚：略薄于颗粒直径；冰太厚信噪比差，太薄颗粒易变形</li>
</ul>
<h3>2.2 数据采集</h3>
<ul>
<li>300 kV 电镜（Titan Krios）+ 直接电子探测器（Gatan K3 / Falcon 4）</li>
<li><strong>超分辨率模式</strong> + 剂量分配（total dose 约 40–60 e⁻/Å²）</li>
<li>每张照片（movie）通常 30–50 帧，用于运动校正</li>
<li>目标：10,000–20,000 张照片（高分辨率结构通常需要百万级颗粒）</li>
</ul>
<h2>3. 预处理（Preprocessing）</h2>
<h3>3.1 运动校正（Motion Correction）</h3>
<p>电子束引起的样品漂移会模糊图像，需要按帧校正：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 4 中的 motion correction</span>
relion_run_motioncorr \\
  --i Micrographs/ \\
  --o MotionCorr/ \\
  --j 8 \\
  --pipeline_control MotionCorr/job001/
</code></pre>
<p>常用软件：<strong>MotionCor2</strong>（经典）、<strong>RELION 内置</strong>、<strong>cryoSPARC Patch Motion Correction</strong>。</p>
<h3>3.2 CTF 估计</h3>
<p>CTF（Contrast Transfer Function，衬度传递函数）描述了电子显微镜成像的相位翻转与欠焦效应：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># CTFFIND4（经典命令行）</span>
ctffind --input_micrograph MotionCorr/mic1.mrc \\
        --output_diag mic1_diag.mrc \\
        --defU 20000 --defV 20000

<span class="hljs-comment"># Gctf（更快的 GPU 版本）</span>
Gctf --input <span class="hljs-string">'MotionCorr/*.mrc'</span> --apix 1.0
</code></pre>
<p>常用软件：<strong>CTFFIND4</strong>、<strong>Gctf</strong>、<strong>cryoSPARC Patch CTF</strong>。</p>
<p><strong>关键判据</strong>：分辨率截断（Resolution limit）通常取 3–4 Å；剔除 astigmatism 过大或漂移严重的照片。</p>
<h2>4. 颗粒挑选（Particle Picking）</h2>
<h3>4.1 传统方法</h3>
<ul>
<li><strong>模板匹配</strong>（Template-based）：用 2D 类平均做模板，相关性搜索</li>
<li>软件：RELION 的 auto-picking、cryoSPARC 的 Template Picker</li>
</ul>
<h3>4.2 深度学习方法（当前主流）</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Topaz：基于 CNN 的颗粒挑选</span>
topaz train --model-path model.pkl --epochs 20 \\
  --num-workers 4 --train-data particles.toml

topaz extract --model-path model.pkl \\
  --micrographs micrographs.toml \\
  --output-prefix picked --radius 100
</code></pre>
<p>主流工具：<strong>Topaz</strong>、<strong>crYOLO</strong>（速度极快，支持实时挑选）、<strong>cryoSPARC Blob Picker</strong>。</p>
<blockquote>
<p>现代流程常先用 Blob/模板快速挑一批做 2D 分类，再用干净类平均训练 Topaz/crYOLO 模型，显著提升复杂样品（膜蛋白、小蛋白）的挑选质量。</p>
</blockquote>
<h2>5. 2D 分类：筛选颗粒</h2>
<p>2D 分类将颗粒按投影方向聚类，同时<strong>剔除冰污染、错误挑选、变性颗粒</strong>：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 2D classification</span>
relion_refine --o Class2D/ \\
  --particle_images particles.star \\
  --ctf --K 100 --iter 25 --tau2_fudge 4 \\
  --flatten_solvent --zero_mask --pool 30 \\
  --per_particle_ctf --j 8
</code></pre>
<p><strong>质量判据</strong>：</p>
<ul>
<li>类平均应显示清晰的二级结构（α 螺旋、β 片层）</li>
<li>保留占总量 50–80% 的"好"颗粒</li>
<li>好的类平均是后续 3D 重构质量的前提</li>
</ul>
<h2>6. 3D 重构：初始模型与分类</h2>
<h3>6.1 初始模型（Ab-initio）</h3>
<p>无需参考模型，从随机起始重构：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># cryoSPARC</span>
cryosparc_abinitio ... --num_classes 3
</code></pre>
<ul>
<li>常用：<strong>cryoSPARC Ab-initio</strong>（速度快、稳健）、RELION 3D initial model</li>
<li>生成 2–4 个类，选择结构特征最清晰的类作为参考</li>
</ul>
<h3>6.2 3D 分类：处理构象异质性</h3>
<p>生物样品常存在多个构象（开放/闭合、配体结合/未结合）：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># RELION 3D classification</span>
relion_refine --o Class3D/ \\
  --ref initial.mrc --K 4 --iter 40 \\
  --tau2_fudge 8 --particle_diameter 200
</code></pre>
<p><strong>策略</strong>：先粗分（K=4–6）看构象分布，再对感兴趣构象细分（K=2–3），最终每个 3D 类对应一个构象。</p>
<h2>7. 3D 精修（Refinement）</h2>
<h3>7.1 常规精修</h3>
<pre><code class="hljs language-bash">relion_refine --o Refine3D/ \\
  --ref class3d.mrc --particle_diameter 200 \\
  --auto_mask --solvent_mask \\
  --ctf --per_particle_ctf --j 8
</code></pre>
<p>关键参数：</p>
<ul>
<li><strong>C1 对称性</strong>起步，对称分子可指定 C/D/O 对称（如 C3、D7）</li>
<li><strong>溶剂掩膜</strong>（solvent mask）必须使用，否则 FSC 虚高</li>
<li>局部精修（local refinement）处理柔性区域</li>
</ul>
<h3>7.2 贝叶斯抛光（Bayesian Polishing）</h3>
<p>对 movie 帧级颗粒做轨迹拟合与加权，提升高频信息：</p>
<pre><code class="hljs language-bash">relion_polish --i Refine3D/run_data.star \\
  --o Polish/ \\
  --model_type 2 --nr_iter 5
</code></pre>
<h3>7.3 各向异性校正（CTF refinement）</h3>
<p>精修每个颗粒的欠焦、像散、束倾斜：</p>
<pre><code class="hljs language-bash">relion_ctf_refine --i Polish/particles.star \\
  --o CtfRefine/ \\
  --fit_defocus --fit_magnification --fit_tilt
</code></pre>
<h2>8. 分辨率评估：FSC 曲线</h2>
<p><strong>FSC（Fourier Shell Correlation）</strong>：把颗粒随机分成两个半集（half-maps），分别重构后比较两图的相关系数：</p>
<ul>
<li><strong>FSC = 0.143 准则</strong>：金标准分辨率判定（半图校正后）</li>
<li><strong>FSC = 0.5</strong>：保守估计（未校正）</li>
</ul>
<pre><code class="hljs language-bash">relion_postprocess --i Refine3D/run_half1_class001.mrc \\
  --mask mask.mrc --autob_masking \\
  --angpix 1.0
</code></pre>
<p><strong>质量检查表</strong>：</p>
<table>
<thead>
<tr>
<th>指标</th>
<th>良好标准</th>
</tr>
</thead>
<tbody>
<tr>
<td>全局分辨率（FSC 0.143）</td>
<td>≤ 3.5 Å（高分辨）</td>
</tr>
<tr>
<td>颗粒数</td>
<td>≥ 10 万（3 Å 级通常需数十万）</td>
</tr>
<tr>
<td>方向分布</td>
<td>均匀（无 Preferred orientation）</td>
</tr>
<tr>
<td>密度连续性</td>
<td>侧链密度清晰可见（&#x3C; 3 Å）</td>
</tr>
</tbody>
</table>
<h2>9. 常用软件栈对比</h2>
<table>
<thead>
<tr>
<th>环节</th>
<th>RELION（命令行）</th>
<th>cryoSPARC（GUI）</th>
</tr>
</thead>
<tbody>
<tr>
<td>运动校正</td>
<td>内置 / MotionCor2</td>
<td>Patch Motion</td>
</tr>
<tr>
<td>CTF</td>
<td>CTFFIND4 / Gctf</td>
<td>Patch CTF</td>
</tr>
<tr>
<td>挑选</td>
<td>模板 / Topaz 集成</td>
<td>Blob / Template / Topaz</td>
</tr>
<tr>
<td>2D/3D</td>
<td>稳定可靠</td>
<td>快、交互好</td>
</tr>
<tr>
<td>非均匀精修</td>
<td>—</td>
<td><strong>NU-refinement</strong>（柔性大分子利器）</td>
</tr>
</tbody>
</table>
<blockquote>
<p><strong>工作流建议</strong>：cryoSPARC 做预处理与初始重构，RELION 做高精度精修与抛光；或全流程 cryoSPARC（对新手更友好）。</p>
</blockquote>
<h2>10. 小结</h2>
<ul>
<li>预处理（运动校正 + CTF）决定数据质量上限</li>
<li>深度学习颗粒挑选（Topaz/crYOLO）显著提升复杂样品表现</li>
<li>2D 分类筛颗粒，3D 分类处理构象异质性</li>
<li>精修三件套：溶剂掩膜 + 贝叶斯抛光 + CTF refinement</li>
<li>FSC = 0.143 是分辨率金标准</li>
</ul>
<p>后续文章将介绍结构可视化（PyMOL/ChimeraX）与原子建模（Coot/phenix），把密度图变成原子模型。</p>`,F_=`<h1>Atomic Modeling and Refinement</h1>
<p>After obtaining a high-quality density map, the next step is to build and refine an atomic model, ultimately yielding publishable, analyzable coordinate files. This article covers the complete workflow from automatic modeling to quality assessment.</p>
<h2>1. Workflow Overview</h2>
<pre><code>Density map (.mrc)
   ↓
① Initial model (automatic modeling / homologous model / AF prediction)
   ↓
② Model placement into density (docking / rigid body)
   ↓
③ Iterative refinement: Coot manual correction ↔ phenix automatic refinement
   ↓
④ Quality assessment (MolProbity / EMRinger / Q-score)
   ↓
⑤ Validation and publication (PDB deposition)
</code></pre>
<h2>2. Initial Model Sources</h2>
<h3>2.1 Automatic Modeling (Mainstream in the Deep Learning Era)</h3>
<ul>
<li><strong>ModelAngelo</strong>: Automatic modeling based on the AlphaFold architecture, building models directly from density maps</li>
<li><strong>DeepTracer</strong>: Fast fully automatic modeling (GPU)</li>
<li><strong>phenix.map_to_model</strong>: Classic automatic backbone/side-chain building</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ModelAngelo usage example</span>
model_angelo --halfmaps half1.mrc half2.mrc \\
  --output-dir model_angelo_out \\
  --num-workers 8

<span class="hljs-comment"># phenix.map_to_model</span>
phenix.map_to_model map.mrc resolution=3.0 \\
  output_model=auto_model.pdb
</code></pre>
<h3>2.2 Homology Modeling and Molecular Replacement</h3>
<ul>
<li>If a homologous structure exists: use <strong>phenix.dock_in_map</strong> to place it into the density</li>
<li>AlphaFold-predicted models: directly used as initial models (commonly used in cryo-EM)</li>
<li>At low resolution (> 5 Å), place rigid bodies first, then gradually refine</li>
</ul>
<h2>3. Coot: Manual Modeling and Correction</h2>
<p>Coot is the industry standard for manual model correction (GUI tool).</p>
<h3>3.1 Common Operations</h3>
<pre><code>File > Open Coordinates...   # Load model
File > Open Map...           # Load density map (FSC correlation calculated automatically)

# Basic corrections
1. Backbone adjustment: drag Cα (Real Space Refine Zone)
2. Side-chain orientation: Rotamers
3. Missing residues: Add Terminal Residue / Build
4. Delete extra: Delete Residue Range
5. Mutate residues: Mutate &#x26; Autofit
6. Ligand building: Ligand Builder / restraints
</code></pre>
<h3>3.2 Key Concept: Real Space Refine</h3>
<pre><code>Calculate > Real Space Refine Zone
</code></pre>
<ul>
<li>Select a density region, and Coot optimizes the geometry and density fit of that region</li>
<li>Repeat: fix → refine → check → fix again</li>
<li>After each refinement round, return to Coot to check problematic regions (Ramachandran outliers, clashes)</li>
</ul>
<h3>3.3 Ligand Modeling</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Generate ligand restraints</span>
phenix.elbow ligand.smiles    <span class="hljs-comment"># Generate chemical restraints from SMILES</span>
<span class="hljs-comment"># Or</span>
phenix.ligand_idealization ligand.cif
</code></pre>
<p>Use <code>Ligand > Ligand Builder</code> in Coot for manual building, then Real Space Refine.</p>
<h2>4. phenix Automatic Refinement</h2>
<h3>4.1 Basic Refinement</h3>
<pre><code class="hljs language-bash">phenix.real_space_refine model.pdb \\
  map.mrc \\
  resolution=3.0 \\
  output.prefix=refined \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<h3>4.2 Refinement Strategies</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Multiple rounds of iteration: global first, then local</span>
phenix.real_space_refine model.pdb map.mrc \\
  resolution=3.0 \\
  strategy=individual_sites+individual_adp \\
  restraints=ligand.cif \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<p>Commonly used strategy combinations:</p>
<ul>
<li><code>rigid_body</code>: rigid-body refinement first at low resolution</li>
<li><code>individual_sites</code>: per-atom coordinate refinement (high resolution)</li>
<li><code>individual_adp</code>: B-factor refinement</li>
<li><code>torsion</code>: torsion-angle refinement (commonly used for proteins)</li>
<li><code>local_grid</code>: local grid search (flexible regions)</li>
</ul>
<h3>4.3 Automatic Correction Loop (NCS and Water Molecules)</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Find water molecules (at high resolution)</span>
phenix.find_peaks_holes map.mrc map.mrc resolution=3.0

<span class="hljs-comment"># NCS restraints (significantly improves quality when symmetric subunits are present)</span>
phenix.real_space_refine ... strategy=individual_sites+individual_adp ncs_search.enabled=True
</code></pre>
<h2>5. Model Quality Assessment</h2>
<h3>5.1 Geometric Quality: MolProbity</h3>
<pre><code class="hljs language-bash">phenix.molprobity refined.pdb
</code></pre>
<p>Core metrics:</p>
<table>
<thead>
<tr>
<th>Metric</th>
<th>Good standard (3 Å)</th>
<th>Meaning</th>
</tr>
</thead>
<tbody>
<tr>
<td>Ramachandran favored</td>
<td>> 95%</td>
<td>Backbone dihedral angle validity</td>
</tr>
<tr>
<td>Ramachandran outliers</td>
<td>&#x3C; 0.5%</td>
<td>Outlier backbone</td>
</tr>
<tr>
<td>Rotamer outliers</td>
<td>&#x3C; 1%</td>
<td>Abnormal side-chain conformations</td>
</tr>
<tr>
<td>Clashscore</td>
<td>&#x3C; 5</td>
<td>Atomic collisions</td>
</tr>
<tr>
<td>MolProbity score</td>
<td>&#x3C; 2</td>
<td>Composite geometry score</td>
</tr>
</tbody>
</table>
<h3>5.2 Density Fit Quality</h3>
<ul>
<li><strong>EMRinger</strong>: assesses side-chain-to-density fit (> 2 indicates good)</li>
<li><strong>Q-score</strong>: assesses local density quality (1.0 = perfect, 0.5 = acceptable)</li>
<li><strong>CC_mask</strong>: model-density correlation coefficient (> 0.8 good)</li>
<li><strong>map-model FSC</strong>: FSC comparison between the model and the two half-maps</li>
</ul>
<pre><code class="hljs language-bash">phenix.map_model_cc refined.pdb map.mrc
phenix.em_ringer refined.pdb map.mrc
</code></pre>
<h3>5.3 Overfitting Detection</h3>
<p><strong>Free component (CC-free)</strong>: correlation computed independently from the half-map used in refinement. Must be monitored during refinement:</p>
<pre><code>Before refinement: CC_work ≈ CC_free
After refinement: CC_work significantly > CC_free → sign of overfitting
</code></pre>
<p><strong>Standard practice</strong>: use only one half-map for refinement (or full map + independent half-map validation), and use the other half-map to compute CC_free.</p>
<h2>6. Complete Refinement Workflow Example</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. Automatic modeling</span>
model_angelo --halfmaps half1.mrc half2.mrc --output-dir auto

<span class="hljs-comment"># 2. Initial refinement</span>
phenix.real_space_refine auto/model.pdb map.mrc \\
  resolution=3.1 strategy=rigid_body+individual_sites \\
  output.prefix=round1

<span class="hljs-comment"># 3. Coot manual correction (interactive)</span>
coot round1_refined.pdb map.mrc

<span class="hljs-comment"># 4. Second refinement round (add B-factors)</span>
phenix.real_space_refine coot_out.pdb map.mrc \\
  resolution=3.1 strategy=individual_sites+individual_adp \\
  output.prefix=round2

<span class="hljs-comment"># 5. Quality assessment</span>
phenix.molprobity round2_refined.pdb
phenix.map_model_cc round2_refined.pdb map.mrc
phenix.em_ringer round2_refined.pdb map.mrc

<span class="hljs-comment"># 6. Iterate until all metrics pass</span>
</code></pre>
<h2>7. PDB Deposition</h2>
<p>Before publication, prepare:</p>
<ol>
<li>Coordinate file (.pdb / .cif)</li>
<li>Density maps (half-maps + full map, .mrc)</li>
<li>Validation reports (MolProbity, EMRinger, etc.)</li>
<li>FSC curves and mask information</li>
</ol>
<p>Submit via <a href="https://deposit.wwpdb.org/">wwPDB deposition</a> to obtain PDB ID and EMDB ID.</p>
<h2>8. Common Problems</h2>
<table>
<thead>
<tr>
<th>Problem</th>
<th>Cause</th>
<th>Solution</th>
</tr>
</thead>
<tbody>
<tr>
<td>Blurred side-chain density</td>
<td>Insufficient resolution / incorrect B-factors</td>
<td>Increase resolution, check orientation distribution</td>
</tr>
<tr>
<td>Many Ramachandran outliers</td>
<td>Incorrect backbone building</td>
<td>Rebuild the region in Coot based on density</td>
</tr>
<tr>
<td>Overfitting</td>
<td>Excessive refinement</td>
<td>Monitor with CC_free, reduce iterations</td>
</tr>
<tr>
<td>Poor ligand fit</td>
<td>Incomplete restraints</td>
<td>Regenerate with phenix.elbow</td>
</tr>
<tr>
<td>Poor local region quality</td>
<td>Flexibility/conformational heterogeneity</td>
<td>Separate conformations with 3D classification, local refinement</td>
</tr>
</tbody>
</table>
<h2>9. Summary</h2>
<ul>
<li>Three sources for model building: automatic modeling (ModelAngelo/DeepTracer), homology/AF prediction, de novo building</li>
<li>Alternating iterations of Coot manual correction and phenix automatic refinement is the standard workflow</li>
<li>Triple quality check: geometry (MolProbity) + fit (EMRinger/Q-score) + overfitting (CC_free)</li>
<li>Sub-3 Å maps should aim for side-chain-level confidence; ligands and modifications require separate validation</li>
</ul>
<p>This completes the structural biology pipeline: data processing → technical review → visualization → atomic modeling, forming a complete loop from experiment to model.</p>`,$_=`<h1>原子建模与精修</h1>
<p>拿到高质量密度图后，下一步是搭建原子模型并精修，最终获得可发布、可分析的坐标文件。本文覆盖从自动建模到质量评估的完整工作流。</p>
<h2>1. 工作流总览</h2>
<pre><code>密度图（.mrc）
   ↓
① 初始模型（自动建模 / 同源模型 / AF 预测）
   ↓
② 模型放入密度（docking / rigid body）
   ↓
③ 迭代精修：Coot 手动修正 ↔ phenix 自动精修
   ↓
④ 质量评估（MolProbity / EMRinger / Q-score）
   ↓
⑤ 验证与发布（PDB deposition）
</code></pre>
<h2>2. 初始模型来源</h2>
<h3>2.1 自动建模（深度学习时代的主流）</h3>
<ul>
<li><strong>ModelAngelo</strong>：基于 AlphaFold 架构的自动建模，直接从密度图搭模型</li>
<li><strong>DeepTracer</strong>：快速全自动建模（GPU）</li>
<li><strong>phenix.map_to_model</strong>：经典自动建主链/侧链</li>
</ul>
<pre><code class="hljs language-bash"><span class="hljs-comment"># ModelAngelo 使用示例</span>
model_angelo --halfmaps half1.mrc half2.mrc \\
  --output-dir model_angelo_out \\
  --num-workers 8

<span class="hljs-comment"># phenix.map_to_model</span>
phenix.map_to_model map.mrc resolution=3.0 \\
  output_model=auto_model.pdb
</code></pre>
<h3>2.2 同源建模与分子替换</h3>
<ul>
<li>有同源结构：用 <strong>phenix.dock_in_map</strong> 放入密度</li>
<li>AlphaFold 预测模型：直接作为初始模型（cryo-EM 中常用）</li>
<li>低分辨率（> 5 Å）时先放刚性体，再逐步细化</li>
</ul>
<h2>3. Coot：手动建模与修正</h2>
<p>Coot 是手动模型修正的行业标准（GUI 工具）。</p>
<h3>3.1 常用操作</h3>
<pre><code>File > Open Coordinates...   # 载入模型
File > Open Map...           # 载入密度图（自动算 FSC 相关）

# 基本修正
1. 主链调整：拖拽 Cα（Real Space Refine 区域）
2. 侧链方向：Rotamers（旋转异构体）
3. 缺失残基：Add Terminal Residue / Build
4. 删除多余：Delete Residue Range
5. 突变残基：Mutate &#x26; Autofit
6. 配体搭建：Ligand Builder / restraints
</code></pre>
<h3>3.2 关键理念：Real Space Refine</h3>
<pre><code>Calculate > Real Space Refine Zone
</code></pre>
<ul>
<li>选中密度区域，Coot 优化该区域几何与密度拟合</li>
<li>重复：修正 → 精修 → 检查 → 再修正</li>
<li>每轮精修后回到 Coot 检查不良区域（Ramachandran 离群点、clash）</li>
</ul>
<h3>3.3 配体建模</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 生成配体 restraints</span>
phenix.elbow ligand.smiles    <span class="hljs-comment"># 从 SMILES 生成化学约束</span>
<span class="hljs-comment"># 或</span>
phenix.ligand_idealization ligand.cif
</code></pre>
<p>Coot 中 <code>Ligand > Ligand Builder</code> 手动搭建，再 Real Space Refine。</p>
<h2>4. phenix 自动精修</h2>
<h3>4.1 基本精修</h3>
<pre><code class="hljs language-bash">phenix.real_space_refine model.pdb \\
  map.mrc \\
  resolution=3.0 \\
  output.prefix=refined \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<h3>4.2 精修策略</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 多轮迭代：先整体后局部</span>
phenix.real_space_refine model.pdb map.mrc \\
  resolution=3.0 \\
  strategy=individual_sites+individual_adp \\
  restraints=ligand.cif \\
  <span class="hljs-built_in">nproc</span>=8
</code></pre>
<p>常用 strategy 组合：</p>
<ul>
<li><code>rigid_body</code>：低分辨率先做刚性体</li>
<li><code>individual_sites</code>：逐原子坐标精修（高分辨率）</li>
<li><code>individual_adp</code>：B 因子精修</li>
<li><code>torsion</code>：扭转角精修（蛋白质常用）</li>
<li><code>local_grid</code>：局部网格搜索（柔性区域）</li>
</ul>
<h3>4.3 自动修正循环（NCS 与水分子）</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 找水分子（高分辨时）</span>
phenix.find_peaks_holes map.mrc map.mrc resolution=3.0

<span class="hljs-comment"># NCS 约束（有对称亚基时显著提升质量）</span>
phenix.real_space_refine ... strategy=individual_sites+individual_adp ncs_search.enabled=True
</code></pre>
<h2>5. 模型质量评估</h2>
<h3>5.1 几何质量：MolProbity</h3>
<pre><code class="hljs language-bash">phenix.molprobity refined.pdb
</code></pre>
<p>核心指标：</p>
<table>
<thead>
<tr>
<th>指标</th>
<th>良好标准（3 Å）</th>
<th>含义</th>
</tr>
</thead>
<tbody>
<tr>
<td>Ramachandran  favored</td>
<td>> 95%</td>
<td>主链二面角合理性</td>
</tr>
<tr>
<td>Ramachandran outliers</td>
<td>&#x3C; 0.5%</td>
<td>离群主链</td>
</tr>
<tr>
<td>Rotamer outliers</td>
<td>&#x3C; 1%</td>
<td>侧链构象异常</td>
</tr>
<tr>
<td>Clashscore</td>
<td>&#x3C; 5</td>
<td>原子碰撞</td>
</tr>
<tr>
<td>MolProbity score</td>
<td>&#x3C; 2</td>
<td>综合几何分数</td>
</tr>
</tbody>
</table>
<h3>5.2 密度拟合质量</h3>
<ul>
<li><strong>EMRinger</strong>：评估侧链与密度的拟合（> 2 表示良好）</li>
<li><strong>Q-score</strong>：评估局部密度质量（1.0 = 完美，0.5 = 可接受）</li>
<li><strong>CC_mask</strong>：模型-密度相关系数（> 0.8 良好）</li>
<li><strong>map-model FSC</strong>：模型与两张半图的 FSC 对比</li>
</ul>
<pre><code class="hljs language-bash">phenix.map_model_cc refined.pdb map.mrc
phenix.em_ringer refined.pdb map.mrc
</code></pre>
<h3>5.3 过拟合检测</h3>
<p><strong>Free component（CC-free）</strong>：与精修中使用的半图独立计算的相关性。精修时必须监控：</p>
<pre><code>精修前 CC_work ≈ CC_free
精修后 CC_work 明显 > CC_free → 过拟合信号
</code></pre>
<p><strong>规范做法</strong>：精修只用一张半图（或全图 + 独立半图验证），用另一张半图计算 CC_free。</p>
<h2>6. 完整精修流程示例</h2>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 1. 自动建模</span>
model_angelo --halfmaps half1.mrc half2.mrc --output-dir auto

<span class="hljs-comment"># 2. 初始精修</span>
phenix.real_space_refine auto/model.pdb map.mrc \\
  resolution=3.1 strategy=rigid_body+individual_sites \\
  output.prefix=round1

<span class="hljs-comment"># 3. Coot 手动修正（交互）</span>
coot round1_refined.pdb map.mrc

<span class="hljs-comment"># 4. 第二轮精修（加 B 因子）</span>
phenix.real_space_refine coot_out.pdb map.mrc \\
  resolution=3.1 strategy=individual_sites+individual_adp \\
  output.prefix=round2

<span class="hljs-comment"># 5. 质量评估</span>
phenix.molprobity round2_refined.pdb
phenix.map_model_cc round2_refined.pdb map.mrc
phenix.em_ringer round2_refined.pdb map.mrc

<span class="hljs-comment"># 6. 迭代直到所有指标达标</span>
</code></pre>
<h2>7. 提交 PDB</h2>
<p>发布前需准备：</p>
<ol>
<li>坐标文件（.pdb / .cif）</li>
<li>密度图（half-maps + full map，.mrc）</li>
<li>验证报告（MolProbity、EMRinger 等）</li>
<li>FSC 曲线与 mask 信息</li>
</ol>
<p>通过 <a href="https://deposit.wwpdb.org/">wwPDB deposition</a> 提交，获取 PDB ID 与 EMDB ID。</p>
<h2>8. 常见问题</h2>
<table>
<thead>
<tr>
<th>问题</th>
<th>原因</th>
<th>解决</th>
</tr>
</thead>
<tbody>
<tr>
<td>侧链密度模糊</td>
<td>分辨率不足 / B 因子错误</td>
<td>提高分辨率、检查方向分布</td>
</tr>
<tr>
<td>Ramachandran 离群多</td>
<td>主链建错</td>
<td>Coot 中按密度重建该区域</td>
</tr>
<tr>
<td>过拟合</td>
<td>精修过度</td>
<td>用 CC_free 监控、减少迭代</td>
</tr>
<tr>
<td>配体拟合差</td>
<td>restraints 不完整</td>
<td>phenix.elbow 重新生成</td>
</tr>
<tr>
<td>局部区域质量差</td>
<td>柔性/构象异质性</td>
<td>3D 分类分离构象、局部精修</td>
</tr>
</tbody>
</table>
<h2>9. 小结</h2>
<ul>
<li>建模三来源：自动建模（ModelAngelo/DeepTracer）、同源/AF 预测、从头搭建</li>
<li>Coot 手动修正与 phenix 自动精修<strong>交替迭代</strong>是标准节奏</li>
<li>质量三重检查：几何（MolProbity）+ 拟合（EMRinger/Q-score）+ 过拟合（CC_free）</li>
<li>亚 3 Å 应追求侧链级可信度，配体与修饰需单独验证</li>
</ul>
<p>至此结构生物学方向完成：数据处理流程 → 技术综述 → 可视化 → 原子建模，构成从实验到模型的完整闭环。</p>`,B_=`<h1>Biomacromolecular Structure Visualization: PyMOL and ChimeraX in Practice</h1>
<p>Structure visualization is an everyday tool in structural biology. Using two major software packages, PyMOL and UCSF ChimeraX, this article covers everything from PDB data retrieval to publication-quality rendering, encompassing all the core operations needed in daily research.</p>
<h2>1. Data Source: The PDB Database</h2>
<h3>1.1 Search and Download</h3>
<p><a href="https://www.rcsb.org/">wwPDB / RCSB PDB</a> is the global repository for structural data:</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># Command-line download (by PDB ID)</span>
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/7A0A.pdb   <span class="hljs-comment"># Cryo-EM structure example</span>
</code></pre>
<h3>1.2 Core Contents of a PDB File</h3>
<pre><code>HEADER    PLANT SEED PROTEIN             08-JAN-81   1CRN
TITLE     CRAMBIN
ATOM      1  N   THR A   1      11.106  16.700  17.083  1.00 20.14           N
ATOM      2  CA  THR A   1      11.579  17.437  17.166  1.00 20.55           C
...
HELIX    1   1 THR A    1  SER A    8  1                                  8
SHEET    1   A 2 ILE A   9  ARG A  12  0
</code></pre>
<ul>
<li><code>ATOM</code> lines: atomic coordinates (residue, chain, x/y/z, B-factor, element)</li>
<li><code>HELIX</code> / <code>SHEET</code>: secondary structure annotations</li>
<li><code>CONECT</code>: covalent connections</li>
<li>The newer <code>mmCIF</code> format (<code>.cif</code>) has become mainstream and carries more complete information</li>
</ul>
<h2>2. PyMOL Basics</h2>
<h3>2.1 Launch and Load</h3>
<pre><code class="hljs language-bash">pymol 1CRN.pdb          <span class="hljs-comment"># Load directly</span>
<span class="hljs-comment"># Or after launching</span>
<span class="hljs-comment"># File > Open</span>
</code></pre>
<h3>2.2 Basic Display Commands</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Load</span>
fetch 1crn               <span class="hljs-comment"># Fetch from the network</span>
load 1CRN.pdb

<span class="hljs-comment"># Display modes</span>
show cartoon             <span class="hljs-comment"># Cartoon (backbone trace)</span>
show sticks              <span class="hljs-comment"># Sticks (atomic details)</span>
show surface             <span class="hljs-comment"># Surface</span>
hide lines               <span class="hljs-comment"># Hide lines</span>

<span class="hljs-comment"># Coloring</span>
color cyan               <span class="hljs-comment"># Overall color</span>
spectrum count           <span class="hljs-comment"># Rainbow colors (by residue number)</span>
color red, resi <span class="hljs-number">1</span>-<span class="hljs-number">10</span>     <span class="hljs-comment"># Color a specific region</span>
</code></pre>
<h3>2.3 Selections and Operations</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Selection syntax</span>
select helix, ss h                   <span class="hljs-comment"># All helices</span>
select sheet, ss s                   <span class="hljs-comment"># All beta sheets</span>
select active_site, resi <span class="hljs-number">15</span>+<span class="hljs-number">20</span>+<span class="hljs-number">45</span>    <span class="hljs-comment"># Specific residues</span>
select ligand, resn HEM              <span class="hljs-comment"># Ligand (heme)</span>

<span class="hljs-comment"># Object operations</span>
show sticks, active_site
color yellow, active_site
zoom active_site
</code></pre>
<h3>2.4 Measurements</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Distance</span>
distance d1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">20</span>/CA

<span class="hljs-comment"># Angle / dihedral</span>
angle a1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA
dihedral dh1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA, /1CRN//A/<span class="hljs-number">18</span>/CA
</code></pre>
<h3>2.5 Interaction Analysis</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Hydrogen bonds</span>
distance hbond, (resn HEM), (resn HIS)

<span class="hljs-comment"># Contact interface</span>
select interface, chain A within <span class="hljs-number">4.5</span> of chain B
</code></pre>
<h2>3. PyMOL Script Batch Processing</h2>
<p>Write commonly used operations as scripts for reproducible execution:</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># render_view.py</span>
fetch 1crn
hide everything
show cartoon
color spectrum, ss
<span class="hljs-built_in">set</span> cartoon_transparency, <span class="hljs-number">0.2</span>

<span class="hljs-comment"># Zoom into the binding site and render</span>
zoom resi <span class="hljs-number">1</span>-<span class="hljs-number">20</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> ray_opaque_background, <span class="hljs-number">0</span>
ray <span class="hljs-number">1200</span>, <span class="hljs-number">900</span>
png 1crn_view.png, dpi=<span class="hljs-number">300</span>
</code></pre>
<pre><code class="hljs language-bash">pymol -cq render_view.py    <span class="hljs-comment"># -c command-line mode, -q quiet</span>
</code></pre>
<h2>4. UCSF ChimeraX Basics</h2>
<p>ChimeraX has a friendly, modern interface (Qt GUI) and supports batch operations via the <code>open</code> command:</p>
<h3>4.1 Opening Structures</h3>
<pre><code class="hljs language-bash">chimerax 1CRN.pdb
<span class="hljs-comment"># Or via command</span>
open 1crn
open 7a0a
</code></pre>
<h3>4.2 Common Commands</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Display modes</span>
cartoon
stick
surface

<span class="hljs-comment"># Coloring</span>
color bychain          <span class="hljs-comment"># Color by chain</span>
color byhetero         <span class="hljs-comment"># Distinguish ligands/ions</span>
color byattribute bfactor   <span class="hljs-comment"># Color by B-factor</span>

<span class="hljs-comment"># Selections</span>
select :<span class="hljs-number">15</span>-<span class="hljs-number">20</span>          <span class="hljs-comment"># Residues 15-20</span>
select /A               <span class="hljs-comment"># Chain A</span>
select ligand
</code></pre>
<h3>4.3 Cryo-EM Density Map Display (ChimeraX's Strength)</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># Open a density map (.mrc)</span>
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc

<span class="hljs-comment"># Adjust isosurface level</span>
volume level <span class="hljs-number">0.01</span>
volume <span class="hljs-comment">#2 level 0.05</span>

<span class="hljs-comment"># Superpose density map and model</span>
<span class="hljs-built_in">open</span> model.pdb
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc
volume <span class="hljs-comment">#2 level 0.02</span>
color <span class="hljs-comment">#1 cornflowerblue</span>
transparency <span class="hljs-comment">#2 30</span>

<span class="hljs-comment"># Local density inspection</span>
volume zone <span class="hljs-comment">#2 near :45-60 radius 5</span>
</code></pre>
<p><strong>ChimeraX's <code>volume zone</code> is the core tool for checking local density quality</strong>: display density around residues to assess whether side chains are discernible.</p>
<h2>5. Key Points for High-Quality Rendering</h2>
<h3>5.1 Color Principles</h3>
<ul>
<li>Cartoon coloring: by chain (<code>bychain</code>) or by domain</li>
<li>Key residues: use a single highlight color (yellow/red/blue), avoid overusing rainbows</li>
<li>Ligands: <code>byhetero</code> automatically assigns element-based colors (C gray, N blue, O red, S yellow)</li>
</ul>
<h3>5.2 Lighting and Materials</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> specular, <span class="hljs-number">0.5</span>
<span class="hljs-built_in">set</span> ambient, <span class="hljs-number">0.3</span>

<span class="hljs-comment"># ChimeraX</span>
<span class="hljs-built_in">set</span> lightMode full
graphics silhouettes true
</code></pre>
<h3>5.3 Resolution Output</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># PyMOL publication-quality rendering</span>
ray 2400, 1800
png figure.png, dpi=600

<span class="hljs-comment"># ChimeraX</span>
save figure.png width 2400 height 1800 supersample 3
</code></pre>
<h2>6. Structure Comparison and Superposition</h2>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL: align two structures</span>
align model2, model1

<span class="hljs-comment"># ChimeraX: matchmaker</span>
matchmaker <span class="hljs-comment">#2 to #1</span>
</code></pre>
<p><strong>RMSD</strong> is the standard metric for measuring structural similarity (typically reported as Cα or all-atom RMSD).</p>
<h2>7. Interactive Viewing and Sharing</h2>
<ul>
<li><em><em>Mol</em> / NGL</em>*: web-based 3D viewing (built into PDBe, RCSB)</li>
<li><strong>PyMOL Web</strong>: export sessions as interactive HTML pages</li>
<li>Commonly used for scientific communication: PyMOL sessions (<code>.pse</code>), ChimeraX sessions (<code>.cxs</code>)</li>
</ul>
<h2>8. Summary</h2>
<ul>
<li>PDB is the data source; mmCIF is the new standard</li>
<li>PyMOL: powerful commands, scriptable batch processing (<code>pymol -cq script.py</code>)</li>
<li>ChimeraX: modern GUI + density map handling (<code>volume</code>, <code>volume zone</code>) is irreplaceable</li>
<li>Publication-quality rendering: hide the clutter → highlight the key → high-quality lighting → high-resolution export</li>
</ul>
<p>The next article will cover atomic modeling: how to build and refine atomic models from density maps.</p>`,q_=`<h1>生物大分子结构可视化：PyMOL 与 ChimeraX 实战</h1>
<p>结构可视化是结构生物学的日常工具。本文以 PyMOL 与 UCSF ChimeraX 两大主流软件为例，从 PDB 数据获取到出版级渲染，覆盖日常科研所需的全部核心操作。</p>
<h2>1. 数据来源：PDB 数据库</h2>
<h3>1.1 检索与下载</h3>
<p><a href="https://www.rcsb.org/">wwPDB / RCSB PDB</a> 是全球结构数据仓库：</p>
<pre><code class="hljs language-bash"><span class="hljs-comment"># 命令行下载（按 PDB ID）</span>
wget https://files.rcsb.org/download/1CRN.pdb
wget https://files.rcsb.org/download/7A0A.pdb   <span class="hljs-comment"># 冷冻电镜结构示例</span>
</code></pre>
<h3>1.2 PDB 文件的核心内容</h3>
<pre><code>HEADER    PLANT SEED PROTEIN             08-JAN-81   1CRN
TITLE     CRAMBIN
ATOM      1  N   THR A   1      11.106  16.700  17.083  1.00 20.14           N
ATOM      2  CA  THR A   1      11.579  17.437  17.166  1.00 20.55           C
...
HELIX    1   1 THR A    1  SER A    8  1                                  8
SHEET    1   A 2 ILE A   9  ARG A  12  0
</code></pre>
<ul>
<li><code>ATOM</code> 行：原子坐标（残基、链、x/y/z、B 因子、元素）</li>
<li><code>HELIX</code> / <code>SHEET</code>：二级结构注释</li>
<li><code>CONECT</code>：共价连接</li>
<li>新格式 <code>mmCIF</code>（.cif）已成为主流，信息更完整</li>
</ul>
<h2>2. PyMOL 基础</h2>
<h3>2.1 启动与加载</h3>
<pre><code class="hljs language-bash">pymol 1CRN.pdb          <span class="hljs-comment"># 直接加载</span>
<span class="hljs-comment"># 或启动后</span>
<span class="hljs-comment"># File > Open</span>
</code></pre>
<h3>2.2 基础显示命令</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 载入</span>
fetch 1crn               <span class="hljs-comment"># 从网络获取</span>
load 1CRN.pdb

<span class="hljs-comment"># 显示模式</span>
show cartoon             <span class="hljs-comment"># 卡通（主链走向）</span>
show sticks              <span class="hljs-comment"># 棍状（原子细节）</span>
show surface             <span class="hljs-comment"># 表面</span>
hide lines               <span class="hljs-comment"># 隐藏线条</span>

<span class="hljs-comment"># 着色</span>
color cyan               <span class="hljs-comment"># 整体着色</span>
spectrum count           <span class="hljs-comment"># 彩虹色（按残基序号）</span>
color red, resi <span class="hljs-number">1</span>-<span class="hljs-number">10</span>     <span class="hljs-comment"># 指定区域着色</span>
</code></pre>
<h3>2.3 选择与操作</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 选择语法</span>
select helix, ss h                   <span class="hljs-comment"># 所有螺旋</span>
select sheet, ss s                   <span class="hljs-comment"># 所有片层</span>
select active_site, resi <span class="hljs-number">15</span>+<span class="hljs-number">20</span>+<span class="hljs-number">45</span>    <span class="hljs-comment"># 指定残基</span>
select ligand, resn HEM              <span class="hljs-comment"># 配体（血红素）</span>

<span class="hljs-comment"># 对象操作</span>
show sticks, active_site
color yellow, active_site
zoom active_site
</code></pre>
<h3>2.4 测量</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 距离</span>
distance d1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">20</span>/CA

<span class="hljs-comment"># 角度 / 二面角</span>
angle a1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA
dihedral dh1, /1CRN//A/<span class="hljs-number">15</span>/CA, /1CRN//A/<span class="hljs-number">16</span>/CA, /1CRN//A/<span class="hljs-number">17</span>/CA, /1CRN//A/<span class="hljs-number">18</span>/CA
</code></pre>
<h3>2.5 相互作用分析</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 氢键</span>
distance hbond, (resn HEM), (resn HIS)

<span class="hljs-comment"># 接触面</span>
select interface, chain A within <span class="hljs-number">4.5</span> of chain B
</code></pre>
<h2>3. PyMOL 脚本批处理</h2>
<p>把常用操作写成脚本，可重复执行：</p>
<pre><code class="hljs language-python"><span class="hljs-comment"># render_view.py</span>
fetch 1crn
hide everything
show cartoon
color spectrum, ss
<span class="hljs-built_in">set</span> cartoon_transparency, <span class="hljs-number">0.2</span>

<span class="hljs-comment"># 放大结合位点并渲染</span>
zoom resi <span class="hljs-number">1</span>-<span class="hljs-number">20</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> ray_opaque_background, <span class="hljs-number">0</span>
ray <span class="hljs-number">1200</span>, <span class="hljs-number">900</span>
png 1crn_view.png, dpi=<span class="hljs-number">300</span>
</code></pre>
<pre><code class="hljs language-bash">pymol -cq render_view.py    <span class="hljs-comment"># -c 命令行模式，-q 静默</span>
</code></pre>
<h2>4. UCSF ChimeraX 基础</h2>
<p>ChimeraX 界面友好、现代（Qt GUI），且支持 <code>open</code> 命令批量操作：</p>
<h3>4.1 打开结构</h3>
<pre><code class="hljs language-bash">chimerax 1CRN.pdb
<span class="hljs-comment"># 或命令</span>
open 1crn
open 7a0a
</code></pre>
<h3>4.2 常用命令</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 显示模式</span>
cartoon
stick
surface

<span class="hljs-comment"># 配色</span>
color bychain          <span class="hljs-comment"># 按链着色</span>
color byhetero         <span class="hljs-comment"># 配体/离子区分</span>
color byattribute bfactor   <span class="hljs-comment"># 按 B 因子</span>

<span class="hljs-comment"># 选择</span>
select :<span class="hljs-number">15</span>-<span class="hljs-number">20</span>          <span class="hljs-comment"># 残基 15-20</span>
select /A               <span class="hljs-comment"># A 链</span>
select ligand
</code></pre>
<h3>4.3 冷冻电镜密度图展示（ChimeraX 强项）</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># 打开密度图（.mrc）</span>
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc

<span class="hljs-comment"># 调整等值面水平</span>
volume level <span class="hljs-number">0.01</span>
volume <span class="hljs-comment">#2 level 0.05</span>

<span class="hljs-comment"># 密度图与模型叠加</span>
<span class="hljs-built_in">open</span> model.pdb
<span class="hljs-built_in">open</span> <span class="hljs-built_in">map</span>.mrc
volume <span class="hljs-comment">#2 level 0.02</span>
color <span class="hljs-comment">#1 cornflowerblue</span>
transparency <span class="hljs-comment">#2 30</span>

<span class="hljs-comment"># 区域密度查看（局部）</span>
volume zone <span class="hljs-comment">#2 near :45-60 radius 5</span>
</code></pre>
<p><strong>ChimeraX 的 <code>volume zone</code> 是检查局部密度质量的核心工具</strong>：在残基周围显示密度，判断侧链是否可辨。</p>
<h2>5. 高质量渲染要点</h2>
<h3>5.1 色彩原则</h3>
<ul>
<li>卡通着色：按链（bychain）或按结构域</li>
<li>关键残基：统一高亮色（黄/红/蓝），避免彩虹滥用</li>
<li>配体：<code>byhetero</code> 自动区分元素色（C 灰、N 蓝、O 红、S 黄）</li>
</ul>
<h3>5.2 光照与材质</h3>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL</span>
<span class="hljs-built_in">set</span> ray_shadows, <span class="hljs-number">1</span>
<span class="hljs-built_in">set</span> specular, <span class="hljs-number">0.5</span>
<span class="hljs-built_in">set</span> ambient, <span class="hljs-number">0.3</span>

<span class="hljs-comment"># ChimeraX</span>
<span class="hljs-built_in">set</span> lightMode full
graphics silhouettes true
</code></pre>
<h3>5.3 分辨率输出</h3>
<pre><code class="hljs language-bash"><span class="hljs-comment"># PyMOL 出版级渲染</span>
ray 2400, 1800
png figure.png, dpi=600

<span class="hljs-comment"># ChimeraX</span>
save figure.png width 2400 height 1800 supersample 3
</code></pre>
<h2>6. 结构比对与叠合</h2>
<pre><code class="hljs language-python"><span class="hljs-comment"># PyMOL：对齐两个结构</span>
align model2, model1

<span class="hljs-comment"># ChimeraX：matchmaker</span>
matchmaker <span class="hljs-comment">#2 to #1</span>
</code></pre>
<p><strong>RMSD</strong> 是衡量结构相似度的标准指标（通常报告 Cα 或全原子 RMSD）。</p>
<h2>7. 交互式查看与分享</h2>
<ul>
<li><em><em>Mol</em> / NGL</em>*：网页端 3D 查看（PDBe、RCSB 内置）</li>
<li><strong>PyMOL Web</strong>：把会话导出为 HTML 交互页面</li>
<li>科研交流常用：PyMOL session（.pse）、ChimeraX session（.cxs）</li>
</ul>
<h2>8. 小结</h2>
<ul>
<li>PDB 是数据源头，mmCIF 是新标准</li>
<li>PyMOL：命令强大、脚本化批处理（<code>pymol -cq script.py</code>）</li>
<li>ChimeraX：现代 GUI + 密度图处理（<code>volume</code>、<code>volume zone</code>）无可替代</li>
<li>出版级渲染：隐藏杂项 → 突出关键 → 高质量光线 → 高分辨率导出</li>
</ul>
<p>下一篇将介绍原子建模：如何从密度图搭建并精修原子模型。</p>`,z_=(s,n=".json")=>{const e=ee()==="en-US";return s.endsWith(Nc)?`${s}${n}`:e?`${s}${Nc}${n}`:`${s}${n}`},tt=Object.assign({"/data-branch/content/about-en.json":Nb,"/data-branch/content/about.json":$b,"/data-branch/content/categories-en.json":qb,"/data-branch/content/categories.json":Ub,"/data-branch/content/notes-en.json":Vb,"/data-branch/content/notes.json":Wb,"/data-branch/content/posts-en.json":Xb,"/data-branch/content/posts.json":Qb,"/data-branch/content/projects-en.json":Zb,"/data-branch/content/projects.json":n_,"/data-branch/content/resources-en.json":e_,"/data-branch/content/resources.json":l_,"/data-branch/content/tags-en.json":p_,"/data-branch/content/tags.json":r_,"/data-branch/content/topics-en.json":u_,"/data-branch/content/topics.json":d_}),zc=Object.assign({"/data-branch/content/html/notes/ComputationalBiology/protein-design/protein-design-tools-en.html":m_,"/data-branch/content/html/notes/ComputationalBiology/protein-design/protein-design-tools.html":j_,"/data-branch/content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools-en.html":f_,"/data-branch/content/html/notes/ComputationalBiology/virtual-screening/virtual-screening-tools.html":g_,"/data-branch/content/html/notes/Programming/bash/bash-scripting-en.html":b_,"/data-branch/content/html/notes/Programming/bash/bash-scripting.html":__,"/data-branch/content/html/notes/Programming/linux/linux-basics-en.html":y_,"/data-branch/content/html/notes/Programming/linux/linux-basics.html":v_,"/data-branch/content/html/notes/Programming/python/python-advanced-en.html":w_,"/data-branch/content/html/notes/Programming/python/python-advanced.html":C_,"/data-branch/content/html/notes/Programming/python/python-basics-en.html":k_,"/data-branch/content/html/notes/Programming/python/python-basics.html":x_,"/data-branch/content/html/notes/Programming/python/python-data-en.html":P_,"/data-branch/content/html/notes/Programming/python/python-data.html":S_,"/data-branch/content/html/notes/Programming/r/r-basics-en.html":E_,"/data-branch/content/html/notes/Programming/r/r-basics.html":T_,"/data-branch/content/html/notes/Programming/r/r-ggplot2-en.html":A_,"/data-branch/content/html/notes/Programming/r/r-ggplot2.html":R_,"/data-branch/content/html/notes/Programming/r/r-tidyverse-en.html":L_,"/data-branch/content/html/notes/Programming/r/r-tidyverse.html":M_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-overview-en.html":D_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-overview.html":O_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-workflow-en.html":I_,"/data-branch/content/html/notes/StructuralBiology/cryoem/cryoem-workflow.html":N_,"/data-branch/content/html/notes/StructuralBiology/modeling/atomic-modeling-en.html":F_,"/data-branch/content/html/notes/StructuralBiology/modeling/atomic-modeling.html":$_,"/data-branch/content/html/notes/StructuralBiology/visualization/structure-visualization-en.html":B_,"/data-branch/content/html/notes/StructuralBiology/visualization/structure-visualization.html":q_}),Uc=Object.assign({"/data-branch/cache/en/notes/ComputationalBiology/protein-design/protein-design-tools-en.md":()=>Cs(()=>import("./protein-design-tools-en-Ci1PnweF.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/ComputationalBiology/virtual-screening/virtual-screening-tools-en.md":()=>Cs(()=>import("./virtual-screening-tools-en-BN18PK6_.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/bash/bash-scripting-en.md":()=>Cs(()=>import("./bash-scripting-en-DUCfPEss.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/linux/linux-basics-en.md":()=>Cs(()=>import("./linux-basics-en-CmF2DSwL.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-advanced-en.md":()=>Cs(()=>import("./python-advanced-en-D8iy_BGO.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-basics-en.md":()=>Cs(()=>import("./python-basics-en-BREIlyQp.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/python/python-data-en.md":()=>Cs(()=>import("./python-data-en-BL-OlVb4.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-basics-en.md":()=>Cs(()=>import("./r-basics-en-Dmq0v2IN.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-ggplot2-en.md":()=>Cs(()=>import("./r-ggplot2-en-DcjWyXfd.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/Programming/r/r-tidyverse-en.md":()=>Cs(()=>import("./r-tidyverse-en-DUrsIHdv.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/cryoem/cryoem-overview-en.md":()=>Cs(()=>import("./cryoem-overview-en-BjmylUcX.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/cryoem/cryoem-workflow-en.md":()=>Cs(()=>import("./cryoem-workflow-en-Tc1AFAMz.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/modeling/atomic-modeling-en.md":()=>Cs(()=>import("./atomic-modeling-en-DSC3cz6t.js"),[]).then(s=>s.default),"/data-branch/cache/en/notes/StructuralBiology/visualization/structure-visualization-en.md":()=>Cs(()=>import("./structure-visualization-en-BM56MuMO.js"),[]).then(s=>s.default),"/data-branch/content-src/assets/README.md":()=>Cs(()=>import("./README-BHV79cHn.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/protein-design/protein-design-tools.md":()=>Cs(()=>import("./protein-design-tools-C2wDyUIy.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/ComputationalBiology/virtual-screening/virtual-screening-tools.md":()=>Cs(()=>import("./virtual-screening-tools-BdHWOXdu.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/bash/bash-scripting.md":()=>Cs(()=>import("./bash-scripting-DEkzGDwa.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/linux/linux-basics.md":()=>Cs(()=>import("./linux-basics-UyN4mB_u.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-advanced.md":()=>Cs(()=>import("./python-advanced-iN1UORRl.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-basics.md":()=>Cs(()=>import("./python-basics-uh9Cz9FQ.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/python/python-data.md":()=>Cs(()=>import("./python-data-CmFq1YcA.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-basics.md":()=>Cs(()=>import("./r-basics-CqBb8X8p.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-ggplot2.md":()=>Cs(()=>import("./r-ggplot2-CK3q5J49.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/Programming/r/r-tidyverse.md":()=>Cs(()=>import("./r-tidyverse-v7eXFCz-.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/cryoem/cryoem-overview.md":()=>Cs(()=>import("./cryoem-overview-DT2Crx5Y.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/cryoem/cryoem-workflow.md":()=>Cs(()=>import("./cryoem-workflow-C2zs9nW6.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/modeling/atomic-modeling.md":()=>Cs(()=>import("./atomic-modeling-BGRUqNS4.js"),[]).then(s=>s.default),"/data-branch/content-src/notes/StructuralBiology/visualization/structure-visualization.md":()=>Cs(()=>import("./structure-visualization-B3-STJgI.js"),[]).then(s=>s.default),"/data-branch/content-src/templates/new-post.md":()=>Cs(()=>import("./new-post-kZswKwDU.js"),[]).then(s=>s.default)}),te=(s,n)=>{const a=z_(s,".json"),e=Object.keys(tt).find(t=>t.includes(a));if(e)return tt[e].default;if(console.error(`Failed to load JSON content: ${a}`),ee()==="en-US"){const t=Object.keys(tt).find(l=>l.includes(`${s}.json`));if(t)return tt[t].default}return n},U_=s=>{const a=ee()==="en-US",e=s.replace(/\.md$/i,""),t=[];e.endsWith("-en")?(t.push(e),a&&t.push(e.replace(/-en$/,""))):(t.push(a?`${e}-en`:e),t.push(e));for(const l of t){const o=Object.keys(zc).find(p=>p.includes(`/content/html/${l}.html`));if(o)return zc[o]}return`<h1>Content Not Available</h1>
<p>The requested content could not be loaded.</p>`},Yt=()=>te("categories",[]),H_=()=>te("posts",[]),Eu=()=>te("notes",[]),V_=()=>te("tags",[]),G_=()=>te("resources",[]),Tu={introduction:"",experience:[],section:[],contacts:[]},W_=()=>te("about",Tu),K_=async s=>{const n=ee()==="en-US"?"-en":"",a=`${s.replace(/\.md$/,"")}${n}.md`,e=Object.keys(Uc).find(t=>t.endsWith(`/${a}`)||t.endsWith(a));if(!e)return null;try{return await Uc[e]()}catch(t){return console.error(`Failed to load markdown source: ${a}`,t),null}},X_={class:"article-item"},Y_=["aria-expanded"],Q_={class:"tree-item__text"},J_={class:"article-list tree-sublist"},Z_=Es({__name:"TreeNode",props:{dir:{},depth:{},collapsedIds:{},isActive:{type:Function},toArticle:{type:Function}},emits:["toggle"],setup(s,{emit:n}){const a=s,e=n,t=J(()=>!a.collapsedIds.includes(a.dir.id));return(l,o)=>{const p=aa("router-link"),c=aa("TreeNode",!0);return O(),$("li",X_,[x("button",{type:"button",class:Fs(["tree-item tree-item--folder",[s.depth>=2?`tree-item--l${s.depth}`:"",{"tree-item--branch-open":t.value}]]),"aria-expanded":t.value,onClick:o[0]||(o[0]=r=>e("toggle",s.dir))},[x("span",Q_,V(s.dir.name),1),x("i",{class:Fs(["fas fa-chevron-right tree-caret",{"tree-caret--open":t.value}]),"aria-hidden":"true"},null,2)],10,Y_),gn(x("ul",J_,[(O(!0),$(js,null,Is(s.dir.files,r=>(O(),$("li",{key:r.path,class:"article-item"},[fs(p,{to:s.toArticle(r.path),class:Fs(["tree-item",[`tree-item--l${s.depth+1}`,{"tree-item--active":s.isActive(r.path)}]])},{default:an(()=>[ks(V(r.title),1)]),_:2},1032,["to","class"])]))),128)),(O(!0),$(js,null,Is(s.dir.children,r=>(O(),rn(c,{key:r.id,dir:r,depth:s.depth+1,"collapsed-ids":s.collapsedIds,"is-active":s.isActive,"to-article":s.toArticle,onToggle:o[1]||(o[1]=i=>e("toggle",i))},null,8,["dir","depth","collapsed-ids","is-active","to-article"]))),128))],512),[[ii,t.value]])])}}}),sy={class:"navigation-tree"},ny={class:"tree-label"},ay={class:"article-list"},ey=Es({__name:"NavigationTree",setup(s){const n=ra(),a=ps([]),e=ps(""),t=ps([]),{data:l}=wa(()=>Yt(),[]);function o(u){return Tt(e.value)===Tt(u)}function p(u){return!t.value.includes(u.id)}function c(u){t.value=p(u)?[...t.value,u.id]:t.value.filter(h=>h!==u.id)}function r(){const u=e.value;if(!u){a.value=[];return}const h=u.split("/").filter(Boolean);if(h.length<2){a.value=[];return}const d=h[0],f=h[1];let m=null;s:for(const w of l.value)for(const R of w.items){if(R.name!==f)continue;const N=A=>A.articleUrl.includes(`${ta}/${d}/`);if(R.articles?.some(N)||R.categories.some(A=>A.articles.some(T=>T.articleUrl.includes(`${ta}/${d}/${f}/`)))){m=R;break s}}if(!m){a.value=[];return}const k=[],j=[],v=(w,R)=>({title:w,path:eo(R)});m.articles?.forEach(w=>{w.articleUrl&&k.push(v(w.title,w.articleUrl))}),m.categories.forEach(w=>{const R=[];w.articles.forEach(N=>{N.articleUrl&&R.push(v(N.title,N.articleUrl))}),R.length&&j.push({id:`${m.name}/${w.key||w.title}`,name:w.title||w.key,files:R})});const y=[{name:m.title||m.name||f,files:k,children:j}];a.value=y}function i(){e.value=dt(n.params.path),r()}return Ws(()=>n.params.path,i,{immediate:!0}),Ws(l,()=>r()),kn(()=>{r()}),(u,h)=>{const d=aa("router-link");return O(),$("div",sy,[(O(!0),$(js,null,Is(a.value,f=>(O(),$("div",{key:f.name,class:"tree-group"},[x("div",ny,V(f.name),1),x("ul",ay,[(O(!0),$(js,null,Is(f.children,m=>(O(),rn(Z_,{key:m.id,dir:m,depth:1,"collapsed-ids":t.value,"is-active":o,"to-article":Y(Fe),onToggle:c},null,8,["dir","collapsed-ids","to-article"]))),128)),(O(!0),$(js,null,Is(f.files,m=>(O(),$("li",{key:m.path,class:"article-item"},[fs(d,{to:Y(Fe)(m.path),class:Fs(["tree-item tree-item--l1",{"tree-item--active":o(m.path)}])},{default:an(()=>[ks(V(m.title),1)]),_:2},1032,["to","class"])]))),128))])]))),128))])}}}),Qs=(s,n)=>{const a=s.__vccOpts||s;for(const[e,t]of n)a[e]=t;return a},Au=Qs(ey,[["__scopeId","data-v-5d73a6b9"]]),At={search:'<circle cx="11" cy="11" r="8" /><line x1="21" y1="21" x2="16.65" y2="16.65" />',sun:'<circle cx="12" cy="12" r="5" /><line x1="12" y1="1" x2="12" y2="3" /><line x1="12" y1="21" x2="12" y2="23" /><line x1="4.22" y1="4.22" x2="5.64" y2="5.64" /><line x1="18.36" y1="18.36" x2="19.78" y2="19.78" /><line x1="1" y1="12" x2="3" y2="12" /><line x1="21" y1="12" x2="23" y2="12" /><line x1="4.22" y1="19.78" x2="5.64" y2="18.36" /><line x1="18.36" y1="5.64" x2="19.78" y2="4.22" />',moon:'<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z" />',language:'<path d="m5 8 6 6" /><path d="m4 14 6-6 2-3" /><path d="M2 5h12" /><path d="M7 2h1" /><path d="m22 22-5-10-5 10" /><path d="M14 18h6" />',menu:'<line x1="3" y1="12" x2="21" y2="12" /><line x1="3" y1="6" x2="21" y2="6" /><line x1="3" y1="18" x2="21" y2="18" />',github:'<path d="M9 19c-5 1.5-5-2.5-7-3m14 6v-3.87a3.37 3.37 0 0 0-.94-2.61c3.14-.35 6.44-1.54 6.44-7A5.44 5.44 0 0 0 20 4.77 5.07 5.07 0 0 0 19.91 1S18.73.65 16 2.48a13.38 13.38 0 0 0-7 0C6.27.65 5.09 1 5.09 1A5.07 5.07 0 0 0 5 4.77a5.44 5.44 0 0 0-1.5 3.78c0 5.42 3.3 6.61 6.44 7A3.37 3.37 0 0 0 9 18.13V22" />',mail:'<path d="M4 4h16c1.1 0 2 .9 2 2v12c0 1.1-.9 2-2 2H4c-1.1 0-2-.9-2-2V6c0-1.1.9-2 2-2z" /><polyline points="22,6 12,13 2,6" />',copy:'<rect x="9" y="9" width="13" height="13" rx="2" /><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1" />',check:'<polyline points="20 6 9 17 4 12" />',close:'<line x1="18" y1="6" x2="6" y2="18" /><line x1="6" y1="6" x2="18" y2="18" />',"arrow-up":'<line x1="12" y1="19" x2="12" y2="5" /><polyline points="5 12 12 5 19 12" />',list:'<line x1="8" y1="6" x2="21" y2="6" /><line x1="8" y1="12" x2="21" y2="12" /><line x1="8" y1="18" x2="21" y2="18" /><line x1="3" y1="6" x2="3.01" y2="6" /><line x1="3" y1="12" x2="3.01" y2="12" /><line x1="3" y1="18" x2="3.01" y2="18" />'};function yl(s,n=18){return`<svg width="${n}" height="${n}" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.75" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">${At[s]}</svg>`}const ty=["aria-label"],ly=["width","height","innerHTML"],oy=Es({__name:"IconButton",props:{icon:{},ariaLabel:{},size:{default:18}},setup(s){const n=a=>{a.currentTarget?.blur()};return(a,e)=>(O(),$("button",{type:"button",class:"icon-btn","aria-label":s.ariaLabel,onFocus:n},[(O(),$("svg",{class:"icon-btn__svg",width:s.size,height:s.size,viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true",innerHTML:Y(At)[s.icon]},null,8,ly))],40,ty))}}),lt=Qs(oy,[["__scopeId","data-v-b9ed4622"]]),py=Es({__name:"NavActions",props:{mobile:{type:Boolean,default:!1}},emits:["open-search","toggle-menu"],setup(s,{emit:n}){const a=n,{t:e}=Hs(),t=ra(),l=Ge(),o=Ho(),p=J(()=>o.theme==="dark"?!0:o.theme==="light"?!1:typeof window<"u"&&window.matchMedia("(prefers-color-scheme: dark)").matches),c=()=>o.toggleTheme(),r=()=>Ob(l,t);return(i,u)=>(O(),$("div",{class:Fs(["d-flex",s.mobile?"d-lg-none ms-auto app-nav__actions app-nav__actions--mobile":"d-none d-lg-flex ms-auto app-nav__actions"])},[fs(lt,{icon:"search",ariaLabel:Y(e)("search"),onClick:u[0]||(u[0]=h=>a("open-search"))},null,8,["ariaLabel"]),fs(lt,{icon:p.value?"sun":"moon",ariaLabel:Y(e)("theme"),onClick:c},null,8,["icon","ariaLabel"]),fs(lt,{icon:"language",ariaLabel:Y(e)("language"),onClick:r},null,8,["ariaLabel"]),s.mobile?(O(),rn(lt,{key:0,icon:"menu",ariaLabel:Y(e)("menu"),onClick:u[1]||(u[1]=h=>a("toggle-menu"))},null,8,["ariaLabel"])):is("",!0)],2))}}),Hc=Qs(py,[["__scopeId","data-v-259baced"]]),cy={class:"app-header"},ry={class:"container px-0"},iy={class:"navbar navbar-expand-lg app-nav"},uy={class:"container-fluid d-flex app-nav__inner"},hy={class:"app-nav__wordmark"},dy={class:"navbar-nav mb-2 mb-lg-0 app-nav__links"},my={key:0,class:"app-nav__link-divider","aria-hidden":"true"},jy={class:"offcanvas-panel"},fy={class:"offcanvas-section"},gy={class:"offcanvas-head"},by={class:"offcanvas-card"},_y={class:"list-unstyled m-0"},yy={key:0,class:"offcanvas-section"},vy={class:"offcanvas-head"},wy=Es({__name:"AppHeader",emits:["open-search"],setup(s,{emit:n}){const a=ra(),{t:e}=Hs(),t=Ho(),l=Uo(),o=ps(!1),p=ps(!1),c=n,r=J(()=>[{text:e("categories"),href:ga("/category")},{text:e("resources"),href:ga("/resource")},{text:e("about"),href:ga("/about")}]),i=k=>k===a.path?!0:k!==ga("/")&&a.path.startsWith(k),u=J(()=>a.path.includes(`${ta}/`)),h=()=>{window.innerWidth<992?d():o.value=!o.value},d=()=>{Lb(),p.value=!0},f=()=>{p.value=!1,$c()},m=k=>{const j=k.target;j&&j.closest("a")&&f()};return kn(()=>{t.initTheme(),l.initLocale()}),Dn(()=>{$c()}),(k,j)=>{const v=aa("RouterLink");return O(),$(js,null,[x("header",cy,[x("div",ry,[x("nav",iy,[x("div",uy,[fs(v,{class:"navbar-brand app-nav__brand",to:Y(ga)("/"),onClick:j[0]||(j[0]=y=>o.value=!1)},{default:an(()=>[x("span",hy,[ks(V(Y(jn).author),1),j[4]||(j[4]=x("span",{class:"app-nav__apos"},"’",-1)),j[5]||(j[5]=ks("s blog",-1))])]),_:1},8,["to"]),fs(Hc,{mobile:"",onOpenSearch:j[1]||(j[1]=y=>c("open-search")),onToggleMenu:h}),x("div",{class:Fs(["navbar-collapse collapse",{show:o.value}])},[x("ul",dy,[r.value.length?(O(),$("li",my)):is("",!0),(O(!0),$(js,null,Is(r.value,y=>(O(),$("li",{class:"nav-item",key:y.text},[fs(v,{class:Fs(["nav-link app-nav__link",{active:i(y.href)}]),to:y.href,onClick:j[2]||(j[2]=w=>o.value=!1)},{default:an(()=>[ks(V(y.text),1)]),_:2},1032,["to","class"])]))),128))])],2),fs(Hc,{onOpenSearch:j[3]||(j[3]=y=>c("open-search"))})])])])]),p.value?(O(),$("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:Bn(f,["self"])},[x("div",jy,[x("div",fy,[x("div",gy,V(Y(e)("menu")),1),x("div",by,[x("ul",_y,[(O(!0),$(js,null,Is(r.value,y=>(O(),$("li",{key:y.text,class:"my-1"},[fs(v,{to:y.href,class:Fs(["offcanvas-link d-flex align-items-center",{active:i(y.href)}]),onClick:f},{default:an(()=>[x("span",null,V(y.text),1),j[6]||(j[6]=x("i",{class:"fas fa-chevron-right offcanvas-link__chevron"},null,-1))]),_:2},1032,["to","class"])]))),128))])])]),u.value?(O(),$("div",yy,[x("div",vy,V(Y(e)("tableOfContents")),1),x("div",{class:"offcanvas-tree offcanvas-card",onClick:m},[fs(Au)])])):is("",!0)]),x("div",{class:"offcanvas-backdrop",onClick:f})])):is("",!0)],64)}}}),Cy=Qs(wy,[["__scopeId","data-v-3e6c0786"]]),ky={class:"site-footer"},xy={class:"container"},Py={class:"site-footer__inner"},Sy={class:"footer-copy"},Ey={class:"footer-designed"},Ty={class:"footer-designed__text"},Ay={class:"footer-designed__name"},Ry={class:"footer-designed__icons"},Ly=["href"],My=["innerHTML"],Dy=["href"],Oy=["innerHTML"],Iy=Es({__name:"AppFooter",setup(s){const{t:n}=Hs(),a=new Date().getFullYear(),e=jn.startYear&&jn.startYear<a?`${jn.startYear} - ${a}`:`${a}`;return(t,l)=>(O(),$("footer",ky,[x("div",xy,[x("div",Py,[x("p",Sy,"© "+V(Y(e))+" "+V(Y(jn).author),1),x("div",Ey,[x("p",Ty,[ks(V(Y(n)("designedByPrefix")),1),x("strong",Ay,V(Y(jn).author),1),ks(V(Y(n)("designedBySuffix")),1)]),x("div",Ry,[x("a",{href:Y(jn).github,target:"_blank",rel:"noopener noreferrer",class:"footer-link","aria-label":"GitHub"},[(O(),$("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true",innerHTML:Y(At).github},null,8,My))],8,Ly),x("a",{href:`mailto:${Y(jn).email}`,class:"footer-link","aria-label":"Email"},[(O(),$("svg",{class:"footer-icon",width:"16",height:"16",viewBox:"0 0 24 24",fill:"none",stroke:"currentColor","stroke-width":"1.75","stroke-linecap":"round","stroke-linejoin":"round","aria-hidden":"true",innerHTML:Y(At).mail},null,8,Oy))],8,Dy)])])])])]))}}),Ny=Qs(Iy,[["__scopeId","data-v-c2ed6d91"]]),vl="floating-buttons-base-top";function Fy(s){const{sourceId:n,mode:a,onRelease:e}=s,t=s.buttonHeight??40,l=s.gap??48,o=s.margin??20,p=ps(!1),c=ps(null),r=ps(!1),i=ps(0),u=ps(0),h=ps(s.defaultTop),d=ps(!1);let f=null;const m=()=>({gap:l,minTop:a==="stack"?o:o+l,maxTop:a==="stack"?Math.max(0,window.innerHeight-t-o-l):window.innerHeight-t-o}),k=I=>{const{minTop:D,maxTop:W}=m();return Math.max(D,Math.min(W,I))},j=I=>a==="stack"?I+l:I,v=I=>a==="stack"?I-l:I;function y(){f!==null&&(cancelAnimationFrame(f),f=null),p.value=!1}function w(){c.value=j(h.value),!p.value&&(p.value=!0,f=requestAnimationFrame(()=>{f=null,c.value!==null&&window.dispatchEvent(new CustomEvent(vl,{detail:{baseTop:c.value,source:n}})),p.value=!1}))}function R(I){const D=I.detail,W=D?.baseTop;D?.source!==n&&typeof W=="number"&&(h.value=k(v(W)))}function N(I){I.preventDefault(),r.value=!0,d.value=!1,i.value=I.touches[0].clientY,u.value=h.value}function S(I){d.value=!0;const D=I.touches[0].clientY-i.value;h.value=k(u.value+D),w(),I.preventDefault()}function A(I){I.preventDefault(),r.value=!1,d.value||e()}function T(){window.addEventListener(vl,R)}function z(){y(),window.removeEventListener(vl,R)}return io(y),{isDragging:r,buttonTop:h,clampTop:k,dispatchBaseTop:w,onTouchStart:N,onTouchMove:S,onTouchEnd:A,subscribe:T,unsubscribe:z}}const $y=["aria-label"],By=Es({__name:"FloatingButton",props:{sourceId:{},defaultTop:{},mode:{default:"match"},onRelease:{},ariaLabel:{},show:{type:Boolean,default:!0}},setup(s){const n=s,a=()=>{n.onRelease?.()},{isDragging:e,buttonTop:t,clampTop:l,dispatchBaseTop:o,onTouchStart:p,onTouchMove:c,onTouchEnd:r,subscribe:i,unsubscribe:u}=Fy({sourceId:n.sourceId,defaultTop:n.defaultTop??(typeof window<"u"?window.innerHeight:1024)-100,mode:n.mode,onRelease:a}),h=()=>{e.value||a()};function d(){t.value=l(t.value)}return Ws(()=>n.show,f=>{f&&o()}),kn(()=>{window.addEventListener("resize",d,{passive:!0}),i(),o()}),Dn(()=>{window.removeEventListener("resize",d),u()}),(f,m)=>gn((O(),$("button",{class:"floating-btn d-flex align-items-center justify-content-center",onClick:h,"aria-label":s.ariaLabel,onTouchstart:m[0]||(m[0]=Bn((...k)=>Y(p)&&Y(p)(...k),["prevent","stop"])),onTouchmove:m[1]||(m[1]=Bn((...k)=>Y(c)&&Y(c)(...k),["prevent","stop"])),onTouchend:m[2]||(m[2]=Bn((...k)=>Y(r)&&Y(r)(...k),["prevent","stop"])),style:oa({top:Y(t)+"px"})},[$r(f.$slots,"default",{},void 0)],44,$y)),[[ii,s.show]])}}),Ru=Qs(By,[["__scopeId","data-v-c4c6a9da"]]),qy=Es({__name:"BackToTop",setup(s){const{t:n}=Hs(),a=ps(!1),e=ps(null),t=ps(null),l=(typeof window<"u"?window.innerHeight:1024)-100;function o(){Fa()}function p(){if(e.value)return;const c=document.createElement("div");c.style.cssText="position:absolute;left:0;top:0;width:1px;height:181px;pointer-events:none;visibility:hidden",document.body.appendChild(c),t.value=c,e.value=new IntersectionObserver(([r])=>{a.value=!r.isIntersecting},{threshold:0}),e.value.observe(c)}return kn(()=>{p()}),Dn(()=>{e.value&&(e.value.disconnect(),e.value=null),t.value&&(t.value.remove(),t.value=null)}),(c,r)=>(O(),rn(Ru,{"source-id":"btt","default-top":l,mode:"match","on-release":o,ariaLabel:Y(n)("backToTop"),show:a.value},{default:an(()=>[...r[0]||(r[0]=[x("i",{class:"fas fa-arrow-up"},null,-1)])]),_:1},8,["ariaLabel","show"]))}}),zy={class:"page-wrap"},Uy={class:"main-content"},Hy=Es({__name:"App",setup(s){const n=ps(!1),a=Lh(()=>Cs(()=>import("./SearchModal-BnGoJJQr.js"),__vite__mapDeps([0,1]))),e=ra(),t=Uo(),{locale:l}=Hs(),o=J(()=>{const u=fu(e.path);return u==="/"||u===""?"":u}),p=J(()=>Kt(Na(l.value)));qt({htmlAttrs:{lang:()=>Et[Na(l.value)]},link:J(()=>{const u=o.value.includes(`${ta}/`),h=u?o.value.replace(/-en$/,""):o.value,d=u&&!o.value.endsWith("-en")?`${o.value}-en`:o.value;return[{rel:"canonical",href:`${jn.url}/${p.value}${o.value}`},{rel:"alternate",hreflang:"zh-CN",href:`${jn.url}/zh${h}`},{rel:"alternate",hreflang:"en",href:`${jn.url}/en${d}`},{rel:"alternate",hreflang:"x-default",href:`${jn.url}/zh${h}`}]}),meta:J(()=>[{property:"og:url",content:`${jn.url}/${p.value}${o.value}`}])});const c=u=>{const h=zo(u);h&&t.setLocale(h)},r=u=>u instanceof HTMLElement?u.tagName==="INPUT"||u.tagName==="TEXTAREA"||u.isContentEditable:!1,i=u=>{u.key==="/"&&!r(u.target)&&(u.preventDefault(),n.value=!0)};return Ws(()=>e.fullPath,c),kn(()=>{c(e.fullPath),window.addEventListener("keydown",i)}),Dn(()=>{window.removeEventListener("keydown",i)}),(u,h)=>{const d=aa("router-view");return O(),$("div",zy,[fs(Cy,{onOpenSearch:h[0]||(h[0]=f=>n.value=!0)}),x("main",Uy,[fs(d,null,{default:an(({Component:f})=>[fs(Fd,{name:"view",mode:"out-in"},{default:an(()=>[(O(),rn(zh(f)))]),_:2},1024)]),_:1})]),fs(Ny),fs(qy),n.value?(O(),rn(Y(a),{key:0,onClose:h[1]||(h[1]=f=>n.value=!1)})):is("",!0)])}}}),Vy=Qs(Hy,[["__scopeId","data-v-41a9b5b1"]]),Gy={class:"empty-state"},Wy={class:"empty-state__text"},Ky=Es({__name:"EmptyState",props:{icon:{default:""},showHome:{type:Boolean,default:!1}},setup(s){const{t:n}=Hs(),a=J(()=>ga("/"));return(e,t)=>{const l=aa("router-link"),o=ae("reveal");return gn((O(),$("div",Gy,[s.icon?(O(),$("i",{key:0,class:Fs(["empty-state__icon",s.icon]),"aria-hidden":"true"},null,2)):is("",!0),x("p",Wy,[$r(e.$slots,"default",{},void 0)]),s.showHome?(O(),rn(l,{key:1,to:a.value,class:"empty-state__home"},{default:an(()=>[t[0]||(t[0]=x("i",{class:"fas fa-house"},null,-1)),ks(V(Y(n)("backHome")),1)]),_:1},8,["to"])):is("",!0)])),[[o]])}}}),Xy=Qs(Ky,[["__scopeId","data-v-6be4e082"]]),Yy=["aria-label"],Qy={class:"pagination__side"},Jy=["aria-label"],Zy={class:"pagination__pages"},sv={key:1,class:"page-ellipsis"},nv=["onClick"],av={key:2,class:"page-ellipsis"},ev={class:"pagination__side pagination__side--right"},tv=["aria-label"],lv=Es({name:"PostPagination",__name:"Pagination",props:{currentPage:{},totalPages:{},middlePages:{},showFirstPage:{type:Boolean},showLastPage:{type:Boolean},showFirstEllipsis:{type:Boolean},showLastEllipsis:{type:Boolean},onGoToPage:{type:Function},onPrev:{type:Function},onNext:{type:Function}},setup(s){const{t:n}=Hs();return(a,e)=>s.totalPages>1?(O(),$("nav",{key:0,class:"pagination","aria-label":Y(n)("pagination")},[x("div",Qy,[s.currentPage>1?(O(),$("button",{key:0,class:"page-btn page-btn--nav",onClick:e[0]||(e[0]=t=>s.onPrev?.()),"aria-label":Y(n)("prevPage")},[...e[4]||(e[4]=[x("i",{class:"fas fa-chevron-left"},null,-1)])],8,Jy)):is("",!0)]),x("div",Zy,[s.showFirstPage?(O(),$("button",{key:0,class:Fs(["page-btn",{"page-btn--active":s.currentPage===1}]),onClick:e[1]||(e[1]=t=>s.onGoToPage?.(1))}," 1 ",2)):is("",!0),s.showFirstEllipsis?(O(),$("span",sv,"...")):is("",!0),(O(!0),$(js,null,Is(s.middlePages,t=>(O(),$("button",{key:t,class:Fs(["page-btn",{"page-btn--active":s.currentPage===t}]),onClick:l=>s.onGoToPage?.(t)},V(t),11,nv))),128)),s.showLastEllipsis?(O(),$("span",av,"...")):is("",!0),s.showLastPage&&s.totalPages>1?(O(),$("button",{key:3,class:Fs(["page-btn",{"page-btn--active":s.currentPage===s.totalPages}]),onClick:e[2]||(e[2]=t=>s.onGoToPage?.(s.totalPages))},V(s.totalPages),3)):is("",!0)]),x("div",ev,[s.currentPage<s.totalPages?(O(),$("button",{key:0,class:"page-btn page-btn--nav",onClick:e[3]||(e[3]=t=>s.onNext?.()),"aria-label":Y(n)("nextPage")},[...e[5]||(e[5]=[x("i",{class:"fas fa-chevron-right"},null,-1)])],8,tv)):is("",!0)])],8,Yy)):is("",!0)}}),ov=Qs(lv,[["__scopeId","data-v-2ef55601"]]);function Lu(s){return Math.max(1,Math.round(Number(s)/300))}const pv={class:"post-item__index num","aria-hidden":"true"},cv={class:"post-item__meta"},rv={class:"post-item__cat"},iv={class:"post-item__date"},uv={class:"post-item__words"},hv={class:"post-item__reading"},dv={class:"post-item__title"},mv={class:"post-item__preview"},jv={class:"post-item__tags"},fv=["onClick"],gv=Es({__name:"PostCard",props:{post:{},index:{},ordinal:{},categoryLabel:{},articlePath:{}},emits:["tagClick"],setup(s,{emit:n}){const a=n,{t:e}=Hs();function t(l){return e("postReadingTime",{minutes:Lu(l)})}return(l,o)=>{const p=aa("router-link"),c=ae("reveal");return gn((O(),$("article",{class:"post-item",style:oa({"--reveal-delay":Math.min(s.index,5)*40+"ms"})},[x("span",pv,V(String(s.ordinal).padStart(2,"0")),1),fs(p,{to:s.articlePath,class:"post-item__main"},{default:an(()=>[x("div",cv,[x("span",rv,V(s.categoryLabel),1),o[3]||(o[3]=x("span",{class:"divider-v"},null,-1)),x("span",iv,[o[0]||(o[0]=x("i",{class:"fas fa-calendar-alt"},null,-1)),ks(V(s.post.date),1)]),o[4]||(o[4]=x("span",{class:"divider-v"},null,-1)),x("span",uv,[o[1]||(o[1]=x("i",{class:"fas fa-file-lines"},null,-1)),ks(V(s.post.wordCount)+" "+V(Y(e)("wordUnit")),1)]),o[5]||(o[5]=x("span",{class:"divider-v"},null,-1)),x("span",hv,[o[2]||(o[2]=x("i",{class:"fas fa-clock"},null,-1)),ks(V(t(s.post.wordCount)),1)])]),x("h3",dv,V(s.post.title),1),x("p",mv,V(s.post.preview),1)]),_:1},8,["to"]),x("div",jv,[(O(!0),$(js,null,Is(s.post.tags,r=>(O(),$("span",{key:r,class:"post-item__tag",onClick:i=>a("tagClick",r)},V(r),9,fv))),128))]),fs(p,{to:s.articlePath,class:"post-item__arrow","aria-label":s.post.title},{default:an(()=>[...o[6]||(o[6]=[x("i",{class:"fas fa-arrow-right"},null,-1)])]),_:1},8,["to","aria-label"])],4)),[[c]])}}}),bv=Qs(gv,[["__scopeId","data-v-d85b32d3"]]);function Vo(){const{locale:s}=Hs(),n=ra(),a=Ge();function e(l){qc(a,l,{locale:Na(s.value),query:{...n.query,page:"1"}})}function t(l){qc(a,l,{locale:Na(s.value),query:{from:n.fullPath},scroll:!1})}return{goTag:e,goTagFromArticle:t}}function _v(s,n,a=5){if(n<=a)return Array.from({length:n},(p,c)=>c+1);const e=[s];let t=1;for(;e.length<a&&(s-t>=1&&!e.includes(s-t)&&e.push(s-t),!(e.length>=a));)s+t<=n&&!e.includes(s+t)&&e.push(s+t),t++;e.includes(1)||(e.pop(),e.push(1)),e.length<a&&!e.includes(n)?e.push(n):e.length>=a&&!e.includes(n)&&(e.pop(),e.push(n));const l=e.sort((p,c)=>p-c),o=[];for(let p=0;p<l.length;p++)p>0&&l[p]-l[p-1]>1&&o.push("..."),o.push(l[p]);return o[0]!==1&&o.unshift("..."),o}const yv=Es({__name:"PostList",props:{docs:{},perPage:{default:6}},setup(s){const n=s,{t:a,locale:e}=Hs(),t=ra(),l=Ge(),{goTag:o}=Vo(),p=ps(1),c=ps(5),{data:r}=wa(()=>Eu(),[]),{data:i}=wa(()=>Yt(),[]),u=J(()=>Math.max(1,Math.ceil(n.docs.length/n.perPage))),h=J(()=>{const I=(p.value-1)*n.perPage,D=I+n.perPage;return n.docs.slice(I,D)}),d=J(()=>_v(p.value,u.value,c.value)),f=J(()=>d.value.filter(I=>typeof I=="number"&&I!==1&&I!==u.value)),m=J(()=>d.value.includes(1)),k=J(()=>d.value.includes(u.value)),j=J(()=>d.value[0]==="..."||d.value[1]==="..."),v=J(()=>d.value[d.value.length-2]==="..."),y=J(()=>{const I={};return i.value.forEach(D=>D.items.forEach(W=>W.categories.forEach(ls=>{ls.key&&ls.title&&(I[ls.key]=ls.title)}))),I});function w(I){const D=r.value.find(ls=>ls.title===I.title);let W;if(D&&D.relativePath)W=`notes/${D.relativePath}.md`;else{const ls=I.category[1]||"notes",ms=I.title.toLowerCase().replace(/[^a-z0-9]/g,"-");W=`${ls}/${ms}.md`,e.value==="en-US"&&(W=W.replace(".md","-en.md"))}return Fe(W)}function R(I){p.value=I;const D={...t.query,page:String(I)};l.push({path:t.path,query:D}).catch(()=>{}),Fa()}function N(I){I>=1&&I<=u.value&&R(I)}function S(){p.value>1&&R(p.value-1)}function A(){p.value<u.value&&R(p.value+1)}function T(){c.value=window.innerWidth<480?3:5}function z(I){const[D,W]=I;if(!W)return D;const ls=y.value[W]||W;return`${D} / ${ls}`}return Ws(()=>n.docs,()=>{p.value=1}),Ws(()=>t.query.page,I=>{const D=parseInt(String(I)),W=Number.isFinite(D)&&D>=1?Math.min(D,u.value):1;W!==p.value&&(p.value=W,Fa())}),kn(()=>{const I=parseInt(String(t.query.page));p.value=Number.isFinite(I)&&I>=1?Math.min(I,u.value):1,T(),window.addEventListener("resize",T)}),Dn(()=>{window.removeEventListener("resize",T)}),(I,D)=>(O(),$("div",null,[h.value.length?(O(!0),$(js,{key:0},Is(h.value,(W,ls)=>(O(),rn(bv,{key:W.id,post:W,index:ls,ordinal:(p.value-1)*n.perPage+ls+1,"category-label":z(W.category),"article-path":w(W),onTagClick:Y(o)},null,8,["post","index","ordinal","category-label","article-path","onTagClick"]))),128)):(O(),rn(Xy,{key:1,icon:"fas fa-inbox","show-home":""},{default:an(()=>[ks(V(Y(a)("noPosts")),1)]),_:1})),fs(ov,{"current-page":p.value,"total-pages":u.value,"middle-pages":f.value,"show-first-page":m.value,"show-last-page":k.value,"show-first-ellipsis":j.value,"show-last-ellipsis":v.value,"on-go-to-page":N,"on-prev":S,"on-next":A},null,8,["current-page","total-pages","middle-pages","show-first-page","show-last-page","show-first-ellipsis","show-last-ellipsis"])]))}});function le(s){qt({title:J(()=>{const n=typeof s=="string"?s:s?.value;return n?`${n} - ${Ic}`:Ic})})}function vv(s){return s>=1e6?(s/1e6).toFixed(s%1e6?1:0)+"M":s>=1e3?(s/1e3).toFixed(s%1e3?1:0)+"K":String(s)}function wv(s,n=3){const a=new Map;s.forEach(p=>p.tags.forEach(c=>a.set(c,(a.get(c)||0)+1)));const e=Array.from(a.entries()).map(([p,c])=>({name:p,count:c})).sort((p,c)=>p.name.localeCompare(c.name)),t=e.length,l=Math.ceil(t/n),o=new Map([...e].sort((p,c)=>c.count-p.count).map((p,c)=>[p.name,c]));return e.map(p=>{const c=o.get(p.name)??0;return{...p,level:c<l?"lg":c<l*2?"md":"sm"}})}const Cv={class:"page-section home-section"},kv={class:"hero"},xv={class:"hero__main"},Pv={class:"hero__greeting"},Sv={class:"hero__greeting-mark"},Ev={class:"hero__name"},Tv={class:"hero__bio"},Av={class:"hero__stats"},Rv={class:"hero__stat"},Lv={class:"hero__stat-num num"},Mv={class:"hero__stat-label"},Dv={class:"hero__stat"},Ov={class:"hero__stat-num num"},Iv={class:"hero__stat-label"},Nv={class:"hero__stat"},Fv={class:"hero__stat-num num"},$v={class:"hero__stat-label"},Bv={key:0,class:"tags-row"},qv=["onClick"],zv={class:"posts-header"},Uv={key:0,class:"posts-header__tag-title"},Hv={class:"posts-header__actions"},Vv=["aria-label"],Gv=["aria-label"],Wv=Es({name:"HomeView",__name:"Home",setup(s){const{t:n,locale:a}=Hs();le(n("metaHome"));const e=ra(),t=Ge(),{goTag:l}=Vo(),{data:o}=wa(()=>H_(),[]),{data:p}=wa(()=>V_(),[]),c=jn.author,r=J(()=>typeof e.query.tag=="string"?e.query.tag:""),i=J(()=>typeof e.query.from=="string"?e.query.from:""),u=J(()=>{const y=r.value;return y?o.value.filter(w=>w.tags.includes(y)):o.value}),h=J(()=>wv(o.value)),d=J(()=>o.value.length),f=J(()=>p.value.length),m=J(()=>o.value.reduce((y,w)=>y+w.wordCount,0)),k=J(()=>vv(m.value));function j(){const y={...e.query};delete y.tag,delete y.from,y.page="1",t.push({path:e.path,query:y}).catch(()=>{}),Fa()}function v(){i.value&&t.push(i.value).catch(()=>{})}return Ws(a,(y,w)=>{y!==w&&r.value&&j()}),(y,w)=>{const R=ae("reveal");return O(),$("div",Cv,[x("header",kv,[x("div",xv,[x("span",Pv,[x("span",Sv,V(Y(n)("greetingPrefix")),1),ks(V(Y(n)("greeting")),1)]),x("h1",Ev,V(Y(c)),1),x("p",Tv,V(Y(n)("developer")),1)]),x("div",Av,[x("div",Rv,[x("span",Lv,V(d.value),1),x("span",Mv,V(Y(n)("articles")),1)]),x("div",Dv,[x("span",Ov,V(f.value),1),x("span",Iv,V(Y(n)("tags")),1)]),x("div",Nv,[x("span",Fv,V(k.value),1),x("span",$v,V(Y(n)("words")),1)])])]),r.value?is("",!0):gn((O(),$("div",Bv,[(O(!0),$(js,null,Is(h.value,N=>(O(),$("span",{key:N.name,class:Fs(["tag",`tag--${N.level}`]),onClick:S=>Y(l)(N.name)},V(N.name),11,qv))),128))])),[[R]]),gn((O(),$("div",zv,[x("h2",{class:Fs(["posts-header__title",{"posts-header__title--tag":r.value}])},[r.value?(O(),$("span",Uv,"# "+V(r.value),1)):(O(),$(js,{key:1},[ks(V(Y(n)("recentPosts")),1)],64))],2),x("div",Hv,[i.value?(O(),$("button",{key:0,class:"chip-close",onClick:v,"aria-label":Y(n)("backToArticle")},[w[0]||(w[0]=x("i",{class:"fas fa-arrow-left"},null,-1)),ks(V(Y(n)("backToArticle")),1)],8,Vv)):is("",!0),r.value?(O(),$("button",{key:1,class:"chip-close",onClick:j,"aria-label":Y(n)("close")},[w[1]||(w[1]=x("i",{class:"fas fa-times"},null,-1)),ks(V(Y(n)("clearFilter")),1)],8,Gv)):is("",!0)])])),[[R]]),fs(yv,{docs:u.value,perPage:5},null,8,["docs"])])}}}),Kv=Qs(Wv,[["__scopeId","data-v-3f70d3c5"]]),Mu=/^10\.\d{4,9}\//;function Vc(s){return!s||!s.trim()?"":/^https?:\/\//i.test(s)?s:"https://"+s.replace(/^\/+/,"")}function Gc(s){return!s||!s.trim()?"":/^https?:\/\//i.test(s)?s:Mu.test(s)?"https://doi.org/"+s:"https://"+s.replace(/^\/+/,"")}function Xv(s){return s?s.replace(/^https?:\/\//,"").replace(/\/$/,""):""}function Yv(s){return!!s&&(s.includes("doi.org")||Mu.test(s))}const Qv={class:"page-section category-view"},Jv={class:"category-head"},Zv={class:"article-title"},sw={class:"cat-section__header"},nw={class:"cat-section__title"},aw={class:"cat-section__count"},ew={class:"cat-grid"},tw={class:"cat-card__head"},lw={class:"cat-card__name"},ow={class:"cat-card__ext-links"},pw=["href"],cw=["href"],rw={class:"cat-card__desc"},iw={key:0,class:"cat-card__stats"},uw={key:0,class:"cat-stat"},hw={key:1,class:"cat-stat"},dw={key:2,class:"cat-stat"},mw={key:1,class:"cat-card__stats"},jw={key:0,class:"cat-stat"},fw={key:1,class:"cat-stat"},gw={key:2,class:"cat-stat"},bw={key:2,class:"cat-card__tags"},_w={class:"cat-card__links"},yw=["onClick"],vw=Es({name:"CategoryView",__name:"Category",setup(s){const{t:n}=Hs();le(n("metaCategories"));const a=Ge(),{data:e}=wa(()=>Yt(),[]),t=J(()=>n("categories")),l=J(()=>n("notes")),o=J(()=>n("projects")),p=J(()=>n("topics")),c=J(()=>n("seeMore"));function r(f){return f===l.value?"fa-book-open":f===o.value?"fa-folder-open":f===p.value?"fa-flask":"fa-folder"}function i(f){const m=f.items;if(f.title===l.value){const k=m.reduce((j,v)=>j+v.stats.postsCount,0);return n("countPosts",{count:k})}return f.title===o.value?n("countProjects",{count:m.length}):n("countTopics",{count:m.length})}function u(f){return f.stats.latestDate||""}function h(f){return!!(f.url||f.github||f.doi)}function d(f){if(f.root){a.push(ga(f.root)).catch(k=>{k.name!=="NavigationDuplicated"&&!k.toString().includes("Navigation cancelled")&&console.error("Navigation error:",k)});return}const m=["url","github","doi"];for(const k of m){const j=f[k];if(!j)continue;const v=k==="doi"?Gc(j):Vc(j);if(v){window.open(v,"_blank","noopener,noreferrer");return}}}return(f,m)=>{const k=ae("reveal");return O(),$("div",Qv,[x("div",Jv,[gn((O(),$("h1",Zv,[ks(V(t.value),1)])),[[k]])]),(O(!0),$(js,null,Is(Y(e),(j,v)=>(O(),$("div",{key:v,class:"cat-section"},[gn((O(),$("div",sw,[x("h2",nw,[x("i",{class:Fs([["fas",r(j.title)],"cat-section__icon"]),"aria-hidden":"true"},null,2),ks(" "+V(j.title),1)]),x("span",aw,V(i(j)),1)])),[[k]]),x("div",ew,[(O(!0),$(js,null,Is(j.items,(y,w)=>gn((O(),$("div",{key:w,class:"cat-card",style:oa({"--reveal-delay":Math.min(Number(w),5)*40+"ms"})},[x("div",tw,[x("div",lw,V(y.title),1),x("div",ow,[y.github?(O(),$("a",{key:0,href:Y(Vc)(y.github),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"GitHub"},[...m[0]||(m[0]=[x("i",{class:"fab fa-github"},null,-1)])],8,pw)):is("",!0),y.doi?(O(),$("a",{key:1,href:Y(Gc)(y.doi),target:"_blank",rel:"noopener noreferrer",class:"cat-card__ext-link","aria-label":"DOI"},[...m[1]||(m[1]=[x("i",{class:"fas fa-link"},null,-1)])],8,cw)):is("",!0)])]),x("p",rw,V(y.desc),1),j.title===l.value?(O(),$("div",iw,[y.stats?.postsCount?(O(),$("span",uw,[m[2]||(m[2]=x("i",{class:"fas fa-file-lines"},null,-1)),ks(V(Y(n)("countPosts",{count:y.stats.postsCount})),1)])):is("",!0),y.stats?.totalWords?(O(),$("span",hw,[m[3]||(m[3]=x("i",{class:"fas fa-font"},null,-1)),ks(V(Y(n)("countWords",{count:y.stats.totalWords})),1)])):is("",!0),u(y)?(O(),$("span",dw,[m[4]||(m[4]=x("i",{class:"fas fa-clock"},null,-1)),ks(V(u(y)),1)])):is("",!0)])):(O(),$("div",mw,[y.language?(O(),$("span",jw,[m[5]||(m[5]=x("i",{class:"fas fa-code"},null,-1)),ks(V(y.language),1)])):is("",!0),y.year?(O(),$("span",fw,[m[6]||(m[6]=x("i",{class:"fas fa-calendar"},null,-1)),ks(V(y.year),1)])):is("",!0),y.license?(O(),$("span",gw,[m[7]||(m[7]=x("i",{class:"fas fa-scale-balanced"},null,-1)),ks(V(y.license),1)])):is("",!0)])),Array.isArray(y.tags)&&y.tags.length?(O(),$("div",bw,[(O(!0),$(js,null,Is(y.tags,(R,N)=>(O(),$("span",{key:N,class:"cat-card__tag"},V(R),1))),128))])):is("",!0),x("div",_w,[h(y)||y.root?(O(),$("a",{key:0,class:"cat-card__link",onClick:Bn(R=>d(y),["prevent"])},[ks(V(c.value),1),m[8]||(m[8]=x("i",{class:"fas fa-arrow-right"},null,-1))],8,yw)):is("",!0)])],4)),[[k]])),128))])]))),128))])}}}),ww=Qs(vw,[["__scopeId","data-v-7884b90e"]]),Cw={class:"page-section resource-view"},kw={class:"resource-head"},xw={class:"article-title"},Pw={class:"resource-subtitle"},Sw={class:"res-layout"},Ew={class:"res-sidebar"},Tw={class:"res-group__label"},Aw={class:"res-group__count num"},Rw={class:"res-group__items"},Lw=["onClick"],Mw={class:"res-main"},Dw={class:"res-groups"},Ow={class:"res-group-block__title"},Iw={class:"res-grid"},Nw=["href"],Fw={class:"res-card__head"},$w={class:"res-card__name"},Bw={class:"res-card__ext-links"},qw={key:0,class:"res-card__ext-link","aria-label":"DOI"},zw={class:"res-card__desc"},Uw={class:"res-card__footer"},Hw={class:"res-card__url"},Vw=Es({name:"ResourceView",__name:"Resource",setup(s){const{t:n}=Hs();le(n("metaResources"));const{data:a}=wa(()=>G_(),[]),e=ps(null),t=J(()=>n("resources")),l=J(()=>n("resourceSubtitle")),o=J(()=>{const r=e.value;return r?r.children?.length?r.children.filter(i=>i.items?.length):r.items?.length?[{title:r.title,items:r.items}]:[]:[]});function p(r){e.value=r}function c(r){return e.value===r}return Ws(a,r=>{e.value=r[0]?.children?.[0]||null},{immediate:!0}),(r,i)=>{const u=ae("reveal");return O(),$("div",Cw,[x("header",kw,[gn((O(),$("h1",xw,[ks(V(t.value),1)])),[[u]]),x("p",Pw,[i[0]||(i[0]=x("i",{class:"fas fa-circle-info resource-head__icon"},null,-1)),ks(V(l.value),1)])]),x("div",Sw,[x("aside",Ew,[(O(!0),$(js,null,Is(Y(a),h=>(O(),$("div",{key:h.title,class:"res-group"},[x("div",Tw,[x("span",null,V(h.title),1),x("span",Aw,V(h.children?.length||0),1)]),x("div",Rw,[(O(!0),$(js,null,Is(h.children,d=>(O(),$("button",{key:d.title,class:Fs(["res-item",{"res-item--active":c(d)}]),onClick:f=>p(d)},V(d.title),11,Lw))),128))])]))),128))]),i[3]||(i[3]=x("div",{class:"res-divider","aria-hidden":"true"},null,-1)),x("main",Mw,[x("div",Dw,[(O(!0),$(js,null,Is(o.value,h=>(O(),$("section",{key:h.title,class:"res-group-block"},[x("h3",Ow,V(h.title),1),x("div",Iw,[(O(!0),$(js,null,Is(h.items,(d,f)=>gn((O(),$("a",{key:d.name,href:d.url,target:"_blank",rel:"noopener noreferrer",class:"res-card",style:oa({"--reveal-delay":Math.min(Number(f),5)*40+"ms"})},[x("div",Fw,[x("span",$w,V(d.name),1),x("span",Bw,[Y(Yv)(d.url)?(O(),$("span",qw,[...i[1]||(i[1]=[x("i",{class:"fas fa-link"},null,-1)])])):is("",!0)])]),x("p",zw,V(d.desc),1),x("div",Uw,[x("span",Hw,V(Y(Xv)(d.url)),1),i[2]||(i[2]=x("i",{class:"fas fa-arrow-up-right-from-square res-card__arrow"},null,-1))])],12,Nw)),[[u]])),128))])]))),128))])])])])}}}),Gw=Qs(Vw,[["__scopeId","data-v-577a8794"]]),Ww="/assets/avatar-DQvqWlfS.png",Kw={class:"page-section about-view"},Xw={class:"about-head"},Yw={class:"about-head__identity"},Qw={class:"about-head__avatar"},Jw=["src"],Zw={key:1,class:"about-head__initial"},s1={class:"about-head__names"},n1={class:"about-head__name"},a1={class:"about-head__role"},e1={key:0,class:"about-intro"},t1={class:"about-body"},l1={key:0,class:"about-section"},o1={class:"about-section__title"},p1={class:"timeline"},c1={class:"tl-year num"},r1={class:"tl-content"},i1={class:"tl-title"},u1={key:0,class:"tl-desc"},h1={key:1,class:"about-grid"},d1={class:"about-cell__title"},m1={key:0,class:"about-cell__item-name"},j1={key:1,class:"about-cell__item-desc"},f1={class:"about-foot"},g1={class:"about-foot__contacts"},b1=["href"],_1={class:"about-foot__thanks"},y1=Es({name:"AboutView",__name:"About",setup(s){const{t:n}=Hs();le(n("metaAbout"));const{data:a}=wa(()=>W_(),Tu),e=Object.assign({"/src/assets/avatar.png":Ww}),t=J(()=>Object.values(e)[0]||""),l=J(()=>n("thanks")),o=jn.author,p=o.trim().charAt(0).toUpperCase(),c=J(()=>a.value.introduction),r=J(()=>a.value.experience),i=J(()=>a.value.section);return(u,h)=>{const d=ae("reveal");return O(),$("div",Kw,[x("header",Xw,[gn((O(),$("div",Yw,[x("div",Qw,[t.value?(O(),$("img",{key:0,src:t.value,alt:"avatar",class:"about-head__avatar-img",draggable:"false"},null,8,Jw)):(O(),$("span",Zw,V(Y(p)),1))]),x("div",s1,[x("h1",n1,V(Y(o)),1),x("p",a1,V(Y(n)("developer")),1)])])),[[d]]),c.value?gn((O(),$("p",e1,[ks(V(c.value),1)])),[[d]]):is("",!0)]),x("main",t1,[r.value.length?gn((O(),$("section",l1,[x("div",o1,V(Y(n)("experience")),1),x("div",p1,[(O(!0),$(js,null,Is(r.value,(f,m)=>(O(),$("div",{key:m,class:"tl-item"},[x("div",c1,V(f.year),1),x("div",r1,[x("div",i1,V(f.title),1),f.desc?(O(),$("div",u1,V(f.desc),1)):is("",!0)])]))),128))])])),[[d]]):is("",!0),i.value.length?(O(),$("div",h1,[(O(!0),$(js,null,Is(i.value,(f,m)=>gn((O(),$("div",{key:m,class:"about-cell",style:oa({"--reveal-delay":Math.min(Number(m),5)*40+"ms"})},[x("div",d1,V(f.title),1),(O(!0),$(js,null,Is(f.items,(k,j)=>(O(),$("div",{key:j,class:"about-cell__item"},[k.name?(O(),$("span",m1,V(k.name),1)):is("",!0),k.desc?(O(),$("span",j1,V(k.desc),1)):is("",!0)]))),128))],4)),[[d]])),128))])):is("",!0)]),x("footer",f1,[x("div",g1,[(O(!0),$(js,null,Is(Y(a).contacts,f=>(O(),$("a",{key:f.label,href:f.link,target:"_blank",rel:"noopener noreferrer",class:"about-foot__contact"},[x("i",{class:Fs(f.icon)},null,2),x("span",null,V(f.value),1)],8,b1))),128))]),x("p",_1,V(l.value),1)])])}}}),v1=Qs(y1,[["__scopeId","data-v-5ae94f03"]]);function w1(s=1200){const n=ps(!1);let a=null;const e=()=>{a!==null&&(clearTimeout(a),a=null)},t=()=>{e(),n.value=!0,a=setTimeout(()=>{n.value=!1,a=null},s)},l=()=>{e(),n.value=!1};return io(e),{copied:n,showSuccess:t,showFailure:l}}async function Du(s){try{return await navigator.clipboard.writeText(s),!0}catch{const n=document.activeElement,a=document.getSelection(),e=a&&a.rangeCount>0?a.getRangeAt(0).cloneRange():null,t=document.createElement("textarea");t.value=s,t.setAttribute("readonly",""),t.style.position="fixed",t.style.top="0",t.style.left="0",t.style.opacity="0",document.body.appendChild(t),t.focus(),t.select(),t.setSelectionRange(0,t.value.length);let l=!1;try{l=document.execCommand("copy")}catch{l=!1}if(document.body.removeChild(t),n instanceof HTMLElement&&n.focus(),e&&e.startContainer.isConnected){const o=document.getSelection();o?.removeAllRanges(),o?.addRange(e)}return l}}const C1=["innerHTML"],k1=Es({__name:"RenderMarkdown",props:{rawMarkdown:{default:""},articlePath:{default:""},articleTitle:{default:""}},emits:["markdown-rendered"],setup(s,{emit:n}){const{t:a}=Hs(),e=Object.assign({}),t=s,l=n,o=ps(""),p=ct("markdownContainer");async function c(A){const T=r(A,t.articlePath);o.value=T,await vn(),l("markdown-rendered"),d(),f(),i()}function r(A,T){try{const z=T.replace(/^[./]*/,"").replace(/\.md$/,"").split("/").slice(0,-1).join("/"),I=D=>{if(/^(https?:)?\/\//i.test(D)||D.startsWith("/"))return D;const W=(z+"/"+D).split("/").filter(cs=>cs&&cs!=="."),ls=[];W.forEach(cs=>cs===".."?ls.pop():ls.push(cs));const ms=ls.join("/"),us=[`@data/content-src/${ms}`,`${ms}`];for(const cs of us){const $s=Object.keys(e).find(Bs=>Bs.endsWith(`/${ms}`)||Bs===cs);if($s)return e[$s]}const es=D.replace(/^\.\//,"");if(es.startsWith("assets/"))for(const cs of[`@data/content-src/${es}`,es]){const $s=Object.keys(e).find(Bs=>Bs.endsWith(`/${es}`)||Bs===cs);if($s)return e[$s]}return D};return A.replace(/<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi,(D,W,ls,ms)=>`<img ${W}src="${I(ls.trim())}"${ms}>`)}catch(z){return console.warn("rewriteImageLinks failed",z),A}}function i(){const A=p.value;A&&(u(A),h(A))}function u(A){if(!t.articleTitle)return;const T=t.articleTitle.trim().toLowerCase();A.querySelectorAll("h1").forEach(z=>{z.textContent.trim().toLowerCase()===T?z.remove():z.replaceWith(Object.assign(document.createElement("h2"),{...Object.fromEntries(Array.from(z.attributes).map(D=>[D.name,D.value])),innerHTML:z.innerHTML}))})}function h(A){const T=(z,I)=>{_u(z,qo),setTimeout(()=>I.blur(),300)};A.querySelectorAll("h2, h3, h4, h5, h6").forEach(z=>{z.querySelector(".heading-anchor")?.remove();const I=Object.assign(document.createElement("button"),{type:"button",className:"heading-anchor",textContent:"#",ariaLabel:a("anchorHeading"),tabIndex:0,ariaHidden:"false"});I.addEventListener("click",D=>{D.stopPropagation(),T(z,I)}),I.addEventListener("keydown",D=>{(D.key==="Enter"||D.key===" ")&&(D.preventDefault(),T(z,I))}),z.appendChild(I)})}function d(){const A=p.value;A&&A.querySelectorAll("pre").forEach(T=>{if(T.querySelector(".code-block-header"))return;const z=T.querySelector("code");if(!z)return;const I=(z.className.match(/language-(\w+)/)||["","text"])[1],D=document.createElement("div");D.className="code-block-header d-flex align-items-center justify-content-between";const W=document.createElement("span");W.className="code-language",W.textContent=I;const ls=document.createElement("button");ls.type="button",ls.className="copy-button btn-icon d-flex align-items-center justify-content-center",ls.setAttribute("aria-label",a("copyCode")),ls.innerHTML=yl("copy",16),ls.addEventListener("click",()=>R(z.textContent??"",ls)),D.append(W,ls);const ms=document.createElement("div");ms.className="code-block-wrapper",T.parentNode?.insertBefore(ms,T),ms.append(D,T)})}function f(){const A=p.value;A&&A.querySelectorAll("table").forEach(T=>{if(T.closest(".table-copyable"))return;const z=document.createElement("div");z.className="table-copyable";const I=document.createElement("button");I.type="button",I.className="table-copy-btn",I.setAttribute("aria-label",a("copyTable")),I.innerHTML=yl("copy",16),I.addEventListener("click",()=>m(T,I)),z.append(I),T.parentNode?.insertBefore(z,T),z.append(T)})}function m(A,T){const I=Array.from(A.querySelectorAll("tr")).map(D=>Array.from(D.querySelectorAll("th, td")).map(W=>(W.textContent||"").trim().replace(/\s+/g," ")).join("	"));R(I.join(`
`),T)}const{copied:k,showSuccess:j,showFailure:v}=w1(),y=ps(null),w=ps("");async function R(A,T){try{if(!await Du(A)){console.warn(a("copyFailed")),v();return}S(),y.value=T,w.value=T.innerHTML,N(T),j()}catch(z){console.error(a("copyFailed"),z),v()}}function N(A){A.style.color="var(--primary)",A.innerHTML=yl("check",16)}function S(){const A=y.value;A&&(A.innerHTML=w.value,A.style.color="",y.value=null)}return Ws(k,A=>{A||S()}),Ws(()=>t.rawMarkdown,A=>c(A),{immediate:!0}),(A,T)=>(O(),$("div",{class:"markdown-body",innerHTML:o.value,ref_key:"markdownContainer",ref:p},null,8,C1))}}),x1={class:"on-this-page"},P1={class:"otp-header"},S1={class:"otp-title"},E1={class:"otp-list"},T1=["onClick","onKeydown"],A1={class:"otp-text"},R1={key:0,class:"otp-sublist"},L1=["onClick","onKeydown"],M1={class:"otp-subtext"},D1=36,O1=Es({__name:"OnThisPage",props:{containerSelector:{default:".markdown-body"},levels:{default:()=>[2,3,4,5,6]},offset:{default:8}},emits:["navigate"],setup(s,{expose:n,emit:a}){const e=s,t=a,{t:l}=Hs(),o=ps([]),p=ps(""),c=ps(null),r=ps(null),i=ps(null),u=ps(null),h=ps(null),d=ps(null);function f(){i.value&&(i.value.disconnect(),i.value=null),u.value&&(clearTimeout(u.value),u.value=null),h.value&&(clearInterval(h.value),h.value=null),d.value&&(d.value.disconnect(),d.value=null)}function m(){o.value=[],p.value="",f(),vn(()=>v())}function k(){w(),R(),vn(j)}function j(){if(!r.value||!c.value)return;const S=c.value.querySelector(".otp-list");if(!S)return;const A=S.querySelector(".otp-item.active > .otp-link, .otp-subitem.active > .otp-sublink");if(!A){r.value.style.opacity="0";return}r.value.style.top=`${A.offsetTop+D1}px`,r.value.style.opacity="1"}Ws(p,()=>vn(j));function v(){setTimeout(()=>k(),100);const S=()=>{const A=document.querySelector(e.containerSelector);A&&(h.value&&(clearInterval(h.value),h.value=null),i.value||(i.value=new MutationObserver(()=>{clearTimeout(u.value??void 0),u.value=window.setTimeout(()=>k(),100)}),i.value.observe(A,{childList:!0,subtree:!0,attributes:!0,attributeFilter:["id"]})),k())};S(),!h.value&&!document.querySelector(e.containerSelector)&&(h.value=window.setInterval(S,200))}function y(S){try{const A=S.cloneNode(!0);return A.querySelectorAll(".heading-anchor")?.forEach(T=>T.remove()),(A.textContent||"").replace(/\s*#\s*$/,"").trim()}catch{return(S.textContent||"").replace(/\s*#\s*$/,"").trim()}}function w(){const S=document.querySelector(e.containerSelector);if(!S){o.value=[];return}const A=e.levels.map(ms=>`h${ms}`).join(","),T=Array.from(S.querySelectorAll(A));if(T.length===0){o.value=[];return}T.forEach(ms=>{if(!ms.id){let us=ms.textContent.trim().toLowerCase().replace(/[^\u4e00-\u9fa5a-zA-Z0-9\s-]/g,"").replace(/\s+/g,"-");us||(us=`section-${Math.random().toString(36).substring(2,9)}`);let es=us,cs=1;for(;document.getElementById(es);)es=`${us}-${cs++}`;ms.id=es}});const z=new Set(e.levels),I=Math.min(...e.levels),D=I+1,W=[];let ls=null;for(const ms of T){const us=parseInt(ms.tagName.substring(1),10);if(!z.has(us))continue;const es={id:ms.id,text:y(ms),level:us,children:[]};us===I?(W.push(es),ls=es):ls&&us>=D?ls.children.push(es):W.push(es)}o.value=W}function R(){d.value&&(d.value.disconnect(),d.value=null);const S=document.querySelector(e.containerSelector);if(!S)return;const A=e.levels.map(z=>`h${z}`).join(","),T=Array.from(S.querySelectorAll(A));T.length!==0&&(d.value=new IntersectionObserver(z=>{for(const I of z)I.isIntersecting&&(p.value=I.target.id)},{rootMargin:`-${e.offset}px 0px -60% 0px`,threshold:0}),T.forEach(z=>d.value?.observe(z)))}function N(S){t("navigate",S);const A=document.getElementById(S);A&&_u(A,e.offset)}return kn(()=>{w(),R(),vn(()=>{j(),v()})}),Dn(()=>{f()}),n({refreshToc:k,resetToc:m}),(S,A)=>(O(),$("nav",x1,[x("div",{class:"otp-content",ref_key:"otpContent",ref:c},[x("span",{class:"otp-marker",ref_key:"marker",ref:r,"aria-hidden":"true"},null,512),x("div",P1,[x("span",S1,V(Y(l)("tableOfContents")),1)]),x("ul",E1,[(O(!0),$(js,null,Is(o.value,(T,z)=>(O(),$("li",{key:z,class:Fs(["otp-item",{active:p.value===T.id}])},[x("a",{class:"otp-link",role:"button",tabindex:"0",onClick:Bn(I=>N(T.id),["prevent"]),onKeydown:Op(Bn(I=>N(T.id),["prevent"]),["enter"])},[x("span",A1,V(T.text),1)],40,T1),T.children&&T.children.length?(O(),$("ul",R1,[(O(!0),$(js,null,Is(T.children,(I,D)=>(O(),$("li",{key:D,class:Fs(["otp-subitem",{active:p.value===I.id}])},[x("a",{class:"otp-sublink",role:"button",tabindex:"0",onClick:Bn(W=>N(I.id),["prevent"]),onKeydown:Op(Bn(W=>N(I.id),["prevent"]),["enter"])},[x("span",M1,V(I.text),1)],40,L1)],2))),128))])):is("",!0)],2))),128))])],512)]))}}),Ou=Qs(O1,[["__scopeId","data-v-8e56cd8b"]]),I1={class:"offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm"},N1={class:"offcanvas-section"},F1={class:"offcanvas-card"},$1=Es({__name:"TocDrawer",setup(s){const{t:n}=Hs(),a=ps(!1),e=ps(!1),t=ps(null),l=(typeof window<"u"?window.innerHeight:1024)-160;function o(){a.value=window.innerWidth<992}function p(){t.value=Mb(),e.value=!0}function c(){e.value=!1,Bc(t.value)}function r(){c()}return kn(()=>{window.addEventListener("resize",o,{passive:!0}),o()}),Dn(()=>{window.removeEventListener("resize",o),Bc(t.value)}),(i,u)=>(O(),$(js,null,[fs(Ru,{class:"d-lg-none","source-id":"toc","default-top":l,mode:"stack","on-release":p,ariaLabel:Y(n)("openToc"),show:a.value},{default:an(()=>[...u[0]||(u[0]=[x("i",{class:"fas fa-bookmark"},null,-1)])]),_:1},8,["ariaLabel","show"]),e.value?(O(),$("div",{key:0,class:"mobile-offcanvas d-lg-none",onClick:Bn(c,["self"])},[x("div",I1,[x("div",N1,[x("div",F1,[fs(Ou,{containerSelector:".markdown-body",levels:[2,3],offset:Y(qo),onNavigate:r},null,8,["offset"])])])]),x("div",{class:"offcanvas-backdrop",onClick:c})])):is("",!0)],64))}}),B1=Qs($1,[["__scopeId","data-v-2d948131"]]);function q1(s,n){const a=[];return s.forEach(e=>e.items.forEach(t=>{t.articles?.forEach(l=>{a.push(l)}),t.categories.forEach(l=>l.articles.forEach(o=>{a.push(o)}))})),a}const z1={class:"container view-container article-view"},U1={class:"row py-4 px-0"},H1={class:"col-12 col-lg-2 order-2 order-lg-1 docs-sidebar-col"},V1={key:0,class:"navigation-container mb-0"},G1={class:"col-12 col-lg-8 order-1 order-lg-2 docs-main-col",ref:"mainContent"},W1={class:"article-panel"},K1={class:"article-panel__body"},X1={class:"article-content"},Y1={key:0,class:"article-meta"},Q1={class:"article-title"},J1={class:"article-meta__row"},Z1={key:0,class:"meta-line"},s2={key:1,class:"meta-line"},n2=["aria-label"],a2={key:0,class:"article-meta__tags"},e2=["onClick"],t2={key:2,class:"article-navigation"},l2={key:0,class:"article-nav-spacer","aria-hidden":"true"},o2={class:"nav-details"},p2={class:"nav-label"},c2={class:"nav-title"},r2={class:"nav-details"},i2={class:"nav-label"},u2={class:"nav-title"},h2={class:"col-12 col-lg-2 order-3 docs-toc-col"},d2={key:0,class:"toc-container mt-0"},m2=Es({name:"ArticleView",__name:"Article",setup(s){const{t:n,locale:a}=Hs(),e=ra(),{goTagFromArticle:t}=Vo(),l=ps(""),o=ps(""),p=ps([]),c=ps([]),r=ps(typeof window<"u"?window.innerWidth:1024),i=ps(0),u=ps(!1);let h=null,d=!1;const f=ct("onThisPageRef"),m=ct("leftSidebarContent"),k=ct("rightSidebarContent"),j=J(()=>r.value>=992),v=J(()=>!!A.value?.path.startsWith("notes/")),y=J(()=>n("updatedAt")),w=J(()=>n("prevPage")),R=J(()=>n("nextPage")),N=os=>os.replace(/\.md$/,""),S=os=>p.value.find(gs=>Tt(gs.path)===Tt(os))||null,A=J(()=>o.value?S(o.value):null);le(J(()=>A.value?.title)),qt({meta:J(()=>A.value?.description?[{name:"description",content:A.value.description}]:[])});const T=J(()=>{if(!A.value)return[];const[os,gs]=N(A.value.path).split("/"),_s=[],Ds=(Ts,F)=>{if(!F.trim())return;const X=N(eo(F)),[K,as]=X.split("/");K!==os||as!==gs||_s.push({title:Ts,path:`${X}.md`})};for(const Ts of c.value)for(const F of Ts.items)F.name===gs&&(F.articles?.forEach(X=>Ds(X.title,X.articleUrl)),F.categories.forEach(X=>X.articles.forEach(K=>Ds(K.title,K.articleUrl))));return _s}),z=J(()=>{const os=A.value;return os?T.value.findIndex(gs=>N(gs.path)===N(os.path)):-1}),I=J(()=>{const os=z.value;return os>0?T.value[os-1]:null}),D=J(()=>{const os=z.value,gs=T.value.length-1;return os>=0&&os<gs?T.value[os+1]:null}),W=J(()=>Lu(A.value?.wordCount??0));function ls(os){return n("articleReadingTime",{minutes:os})}async function ms(){if(!A.value)return;let os="";try{const _s=await K_(A.value.path);_s&&(os=_s.trim())}catch{}if(!os){const _s=document.querySelector(".markdown-body");if(!_s)return;const Ds=_s.cloneNode(!0);Ds.querySelectorAll(".heading-anchor, .code-block-header, .copy-button, .table-copy-btn").forEach(F=>F.remove());const Ts=(Ds.innerText||"").trim();os=`${A.value.title}

${Ts}`}await Du(os)?(u.value=!0,h&&clearTimeout(h),h=window.setTimeout(()=>{u.value=!1},1200)):console.warn(n("copyFailed"))}function us(){try{const os=Yt();c.value=os;const gs=new Map;os.forEach(_s=>_s.items.forEach(Ds=>{const Ts=Ds.stats.latestDate||"";Ds.articles?.forEach(F=>gs.set(F.articleUrl,Ts)),Ds.categories.forEach(F=>{const X=F.stats.latestDate||Ts;F.articles.forEach(K=>gs.set(K.articleUrl,X))})})),p.value=q1(os).map(_s=>({title:_s.title,path:eo(_s.articleUrl),date:gs.get(_s.articleUrl)??"",tags:_s.tags,wordCount:_s.wordCount}))}catch{p.value=[],c.value=[]}}function es(){r.value=window.innerWidth,Bs()}function cs(){d||(d=!0,requestAnimationFrame(()=>{d=!1;const gs=document.documentElement.scrollHeight-window.innerHeight;i.value=gs>0?Math.min(100,window.scrollY/gs*100):0}))}function $s(){l.value="";try{const os=S(dt(e.params.path));if(!os)throw new Error(`Article not found: ${e.params.path}`);o.value=os.path,l.value=U_(o.value),vn(()=>{typeof window>"u"||(Fa(),Bs(),f.value?.refreshToc())})}catch{l.value=`# Article Not Found

The requested article could not be loaded. Please check the URL.`,vn(()=>{typeof window>"u"||(Fa(),f.value?.refreshToc())})}}function Bs(){if(typeof window>"u")return;const gs=document.querySelector("header")?.offsetHeight||60,_s=window.innerHeight,Ds=Math.max(200,_s-gs-24-24);[m.value,k.value].forEach(F=>{F&&(F.style.maxHeight=`${Ds}px`,F.style.overflowY="auto")})}function qs(){Bs(),vn(()=>f.value?.refreshToc())}return us(),$s(),Ws(a,(os,gs)=>{os!==gs&&(us(),$s())}),Ws(()=>e.params.path,(os,gs)=>{const _s=dt(gs),Ds=dt(os);_s!==Ds&&(f.value?.resetToc(),$s())}),Ws(l,()=>{vn(()=>Bs())}),kn(()=>{r.value=window.innerWidth,window.addEventListener("resize",es),window.addEventListener("scroll",cs,{passive:!0})}),Dn(()=>{window.removeEventListener("resize",es),window.removeEventListener("scroll",cs),h&&clearTimeout(h)}),(os,gs)=>{const _s=aa("router-link");return O(),$("div",z1,[x("div",{class:"reading-progress",style:oa({width:i.value+"%"}),"aria-hidden":"true"},null,4),x("div",U1,[x("div",H1,[x("div",{class:"sticky-sidebar",ref_key:"leftSidebarContent",ref:m},[j.value?(O(),$("div",V1,[fs(Au)])):is("",!0)],512)]),x("div",G1,[x("div",W1,[x("div",K1,[x("div",X1,[A.value?(O(),$("div",Y1,[x("h1",Q1,V(A.value.title),1),x("div",J1,[v.value&&A.value.date?(O(),$("span",Z1,[gs[0]||(gs[0]=x("i",{class:"fas fa-calendar-alt"},null,-1)),ks(V(y.value)+" "+V(A.value.date),1)])):is("",!0),W.value>0?(O(),$("span",s2,[gs[1]||(gs[1]=x("i",{class:"fas fa-clock"},null,-1)),ks(V(ls(W.value)),1)])):is("",!0),x("button",{type:"button",class:Fs(["article-copy-btn",{"article-copy-btn--copied":u.value}]),onClick:ms,"aria-label":Y(n)("copyArticle"),"aria-live":"polite"},[x("i",{class:Fs(u.value?"fas fa-check":"fas fa-copy")},null,2),x("span",null,V(u.value?Y(n)("copied"):Y(n)("copyArticle")),1)],10,n2)]),A.value.tags?.length?(O(),$("div",a2,[(O(!0),$(js,null,Is(A.value.tags,(Ds,Ts)=>(O(),$("span",{key:Ts,class:"article-tag",onClick:F=>Y(t)(Ds)}," # "+V(Ds),9,e2))),128))])):is("",!0)])):is("",!0),l.value?(O(),rn(k1,{key:1,rawMarkdown:l.value,articlePath:A.value?.path||"",articleTitle:A.value?.title||"",onMarkdownRendered:qs},null,8,["rawMarkdown","articlePath","articleTitle"])):is("",!0),l.value?(O(),$("nav",t2,[!I.value&&D.value?(O(),$("span",l2)):is("",!0),I.value?(O(),rn(_s,{key:1,to:Y(Fe)(I.value.path),class:"article-nav-item prev"},{default:an(()=>[gs[2]||(gs[2]=x("div",{class:"nav-arrow"},[x("i",{class:"fas fa-arrow-left"})],-1)),x("div",o2,[x("div",p2,V(w.value),1),x("div",c2,V(I.value.title),1)])]),_:1},8,["to"])):is("",!0),D.value?(O(),rn(_s,{key:2,to:Y(Fe)(D.value.path),class:"article-nav-item next"},{default:an(()=>[x("div",r2,[x("div",i2,V(R.value),1),x("div",u2,V(D.value.title),1)]),gs[3]||(gs[3]=x("div",{class:"nav-arrow"},[x("i",{class:"fas fa-arrow-right"})],-1))]),_:1},8,["to"])):is("",!0)])):is("",!0)])])])],512),x("div",h2,[x("div",{class:"sticky-sidebar",ref_key:"rightSidebarContent",ref:k},[j.value?(O(),$("div",d2,[fs(Ou,{ref_key:"onThisPageRef",ref:f,containerSelector:".markdown-body",levels:[2,3],offset:Y(qo)},null,8,["offset"])])):is("",!0)],512)])]),l.value?(O(),rn(B1,{key:0})):is("",!0)])}}}),j2=Qs(m2,[["__scopeId","data-v-01e8bde1"]]),f2={class:"page-section notfound-view"},g2={class:"notfound-desc"},b2=Es({name:"NotFoundView",__name:"NotFound",setup(s){const{t:n}=Hs();le(n("metaNotFound"));const a=ra(),e=J(()=>`/${a.params.locale==="en"?"en":"zh"}/`);return(t,l)=>{const o=aa("router-link");return O(),$("div",f2,[l[1]||(l[1]=x("h1",{class:"notfound-code"},"404",-1)),x("p",g2,V(Y(n)("pageNotFound")),1),fs(o,{to:e.value,class:"notfound-back"},{default:an(()=>[l[0]||(l[0]=x("i",{class:"fas fa-arrow-left"},null,-1)),ks(V(Y(n)("backHome")),1)]),_:1},8,["to"])])}}}),_2=Qs(b2,[["__scopeId","data-v-fb977b62"]]),y2=s=>[{path:`/${s}/`,name:`${s}-Home`,component:Kv},{path:`/${s}/category`,name:`${s}-Category`,component:ww},{path:`/${s}/resource`,name:`${s}-Resource`,component:Gw},{path:`/${s}/about`,name:`${s}-About`,component:v1},{path:`/${s}${ta}/:path*`,name:`${s}-Article`,component:j2,props:!0,beforeEnter:n=>{const a=n.params.path||[];if(a.length>=3){const e=a[a.length-1],t=a[a.length-2];if(t===e||`${t}-en`===e){const l=[...a.slice(0,-2),e];if(Eu().some(o=>o.relativePath===l.join("/")))return{path:`/${s}${ta}/${l.join("/")}`,replace:!0}}}return!0}}],v2=[{path:"/",redirect:()=>`/${za()}/`},{path:"/zh",redirect:"/zh/"},{path:"/en",redirect:"/en/"},{path:"/category",redirect:s=>`/${za()}${s.path}`},{path:"/resource",redirect:s=>`/${za()}${s.path}`},{path:"/about",redirect:s=>`/${za()}${s.path}`},{path:`${ta}/:path*`,redirect:s=>`/${za()}${s.path}`},{path:"/:pathMatch(.*)*",redirect:s=>`/${za()}${s.path}`},...Pb.map(y2).flat(),{path:"/:locale(zh|en)/:pathMatch(.*)*",name:"NotFound",component:_2}];Cs(()=>import("./bootstrap.esm-inoyDHN5.js"),[]);const w2={mounted(s){if(typeof window>"u"||window.matchMedia("(prefers-reduced-motion: reduce)").matches||!("IntersectionObserver"in window))return;s.classList.add("reveal");const n=new IntersectionObserver(([a])=>{a.isIntersecting&&(s.classList.add("reveal-visible"),n.disconnect())},{threshold:.12});n.observe(s),s.__revealIo=n},unmounted(s){s.__revealIo?.disconnect(),s.__revealIo=null}};Jj(Vy,{routes:v2},({app:s,router:n,initialState:a,isClient:e,routePath:t})=>{const l=hm();s.use(l),a.pinia&&(l.state.value=a.pinia),s.use(n),s.use(Ne),s.mixin(zm),s.directive("reveal",w2);const o=Ho(),p=Uo();o.initTheme(),e&&p.initLocale(),e&&"serviceWorker"in navigator&&window.addEventListener("load",()=>{navigator.serviceWorker.register("/sw.js").catch(()=>{})})});export{ga as A,ta as B,js as F,C2 as T,Cs as _,Dn as a,$c as b,O as c,Es as d,rn as e,x as f,$r as g,Bn as h,Qs as i,Ge as j,an as k,Lb as l,Y as m,vn as n,kn as o,gn as p,zs as q,ps as r,$ as s,V as t,Hs as u,k2 as v,Ws as w,ks as x,is as y,Is as z};
