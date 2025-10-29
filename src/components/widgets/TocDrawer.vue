<template>
  <!-- 悬浮目录按钮（移动端显示） -->
  <button
    v-show="visible"
    class="toc-drawer-btn d-lg-none"
    @click="openDrawer"
    aria-label="打开本页目录"
    @touchstart.prevent.stop="onTouchStart"
    @touchmove.prevent.stop="onTouchMove"
    @touchend.prevent.stop="onTouchEnd"
    :style="{ top: buttonTop + 'px' }"
  >
    <i class="fas fa-bookmark"></i>
  </button>

  <!-- 仅移动端显示的右侧抽屉 -->
  <div v-if="show" class="mobile-offcanvas d-lg-none" @click.self="close">
    <div class="offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm">
      <!-- <div class="offcanvas-header d-flex align-items-center justify-content-end">
        <button class="btn btn-sm p-2 btn-icon" @click="close" @focus="$event.target.blur()" aria-label="关闭">
          <i class="bi bi-x-lg"></i>
        </button>
      </div> -->

      <div class="offcanvas-section">
        <div class="offcanvas-card">
          <OnThisPage
            ref="otp"
            containerSelector=".markdown-body"
            :levels="[2, 3]"
            :offset="8"
          />
        </div>
      </div>
    </div>
    <div class="offcanvas-backdrop" @click="close"></div>
  </div>
</template>

<script>
import OnThisPage from '@/components/layout/OnThisPage.vue';

export default {
  name: 'TocDrawer',
  components: { OnThisPage },
  data() {
    return {
      // 事件源标识与 raf 派发节流
      sourceId: 'toc',
      rafPending: false,
      rafLastBaseTop: null,

      // 按钮状态
      visible: false,
      isDragging: false,
      startY: 0,
      initialTop: 0,
      buttonTop: window.innerHeight - 160,
      touchMoved: false,
      // 抽屉状态
      show: false
    };
  },
  mounted() {
    window.addEventListener('scroll', this.onScroll, { passive: true });
    window.addEventListener('resize', this.onResize, { passive: true });
    window.addEventListener('floating-buttons-base-top', this.syncBaseTop);

    this.onResize();
    this.onScroll();

    // 初始同步：作为“上方按钮”，派发一次基准位置，确保 BackToTop 与之保持固定间距
    const GAP = 48;
    this.rafDispatchBaseTop(this.buttonTop + GAP);
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.onScroll);
    window.removeEventListener('resize', this.onResize);
    window.removeEventListener('floating-buttons-base-top', this.syncBaseTop);
    // 卸载时确保滚动状态恢复
    this.unlockScroll();
  },
  methods: {
    // 工具：统一边界
    getBounds() {
      const BUTTON_HEIGHT = 40;
      const GAP = 48;
      const MARGIN = 20;
      return {
        gap: GAP,
        minTop: MARGIN,
        // 作为上方按钮，需为下方按钮预留 GAP 距离
        maxTop: Math.max(this.$el ? 0 : 0, window.innerHeight - BUTTON_HEIGHT - MARGIN - GAP)
      };
    },
    clampTop(top) {
      const { minTop, maxTop } = this.getBounds();
      return Math.max(minTop, Math.min(maxTop, top));
    },
    // raf 节流派发，避免频繁互相抖动
    rafDispatchBaseTop(baseTop) {
      this.rafLastBaseTop = baseTop;
      if (this.rafPending) return;
      this.rafPending = true;
      requestAnimationFrame(() => {
        window.dispatchEvent(new CustomEvent('floating-buttons-base-top', {
          detail: { baseTop: this.rafLastBaseTop, source: this.sourceId }
        }));
        this.rafPending = false;
      });
    },
    // 接收“下方按钮”的基准位置，使 TOC 保持在其上方固定间距，且不越界
    syncBaseTop(e) {
      const base = e?.detail?.baseTop;
      const source = e?.detail?.source;
      if (source === this.sourceId) return; // 忽略自身派发
      const { gap } = this.getBounds();
      if (typeof base === 'number') {
        const desiredTop = base - gap;
        this.buttonTop = this.clampTop(desiredTop);
      }
    },
    onResize() {
      // 目录按钮在移动端常驻显示
      const isMobile = window.innerWidth < 992;
      this.visible = isMobile;
      // 使用统一边界做收敛，避免越界
      this.buttonTop = this.clampTop(this.buttonTop);
    },
    onScroll() {
      if (!this.isDragging) {
        const isMobile = window.innerWidth < 992;
        this.visible = isMobile;
      }
    },
    openDrawer() {
      if (this.isDragging) return;
      this.lockScroll();
      this.show = true;
    },
    close() {
      this.show = false;
      this.unlockScroll();
    },
    onTouchStart(e) {
      // 阻止默认触摸行为，避免浏览器合成点击导致页面滚动
      e.preventDefault();
      this.isDragging = true;
      this.touchMoved = false;
      this.startY = e.touches[0].clientY;
      this.initialTop = this.buttonTop;
    },
    onTouchMove(e) {
      this.touchMoved = true;
      const currentY = e.touches[0].clientY;
      const diffY = currentY - this.startY;

      let newTop = this.clampTop(this.initialTop + diffY);
      this.buttonTop = newTop;

      // TOC 为上方按钮：派发“下方按钮”的基准位置，让 BackToTop 跟随在 GAP 下方
      const { gap } = this.getBounds();
      this.rafDispatchBaseTop(newTop + gap);

      e.preventDefault();
    },
    onTouchEnd(e) {
      // 阻止默认触摸结束产生的幽灵点击
      e.preventDefault();
      this.isDragging = false;
      if (!this.touchMoved) this.openDrawer();
    },
    lockScroll() {
      // 使用 overflow 锁滚，保留滚动坐标，确保 OnThisPage 导航正常工作
      const docEl = document.documentElement;
      const body = document.body;
      if (docEl) {
        docEl.style.overflow = 'hidden';
        docEl.style.overscrollBehavior = 'contain';
      }
      if (body) {
        body.style.overflow = 'hidden';
        body.style.overscrollBehavior = 'contain';
      }
    },
    unlockScroll() {
      const docEl = document.documentElement;
      const body = document.body;
      if (docEl) {
        docEl.style.overflow = '';
        docEl.style.overscrollBehavior = '';
      }
      if (body) {
        body.style.overflow = '';
        body.style.overscrollBehavior = '';
      }
      // 不主动恢复 scrollTop，避免视觉跳动
    }
  }
};
</script>

<style scoped>
/* 按钮样式 */
.toc-drawer-btn {
  position: fixed;
  right: 30px;
  width: 40px;
  height: 40px;
  background-color: var(--app-btn-bg);
  color: var(--app-text);
  border: 1px solid var(--app-border);
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: var(--app-btn-shadow);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}
.toc-drawer-btn:hover {
  background-color: var(--app-btn-hover-bg);
  box-shadow: var(--app-btn-hover-shadow);
}
.toc-drawer-btn:active {
  transform: scale(0.95);
}

/* 抽屉样式（与 AppHeader 一致），右侧版本 */
.mobile-offcanvas {
  position: fixed;
  inset: 0;
  z-index: 1050;
}
.offcanvas-backdrop {
  position: absolute;
  inset: 0;
  background: var(--app-backdrop-bg);
  z-index: 1;
}
.offcanvas-panel {
  position: absolute;
  top: 0; bottom: 0;
  width: min(85vw, 320px);
  background: var(--app-offcanvas-bg);
  box-shadow: var(--app-offcanvas-shadow);
  overflow-y: auto;
  -webkit-overflow-scrolling: touch;
  padding: 0.75rem 0.75rem 1rem;
  z-index: 2;
}
.offcanvas-right { right: 0; } /* 右侧抽屉 */
.offcanvas-header {
  margin-bottom: 0.5rem;
  padding-bottom: 0.5rem;
}
.offcanvas-section + .offcanvas-section {
  margin-top: 0.75rem;
}
.section-title {
  font-size: 0.95rem;
  color: var(--app-text-muted);
  margin-bottom: 0.25rem;
  font-weight: 600;
}
.offcanvas-card {
  background-color: var(--app-card-bg);
  border: none;
  border-radius: 0.5rem;
  padding: 0.5rem;
  margin: 1rem 0;
}
</style>