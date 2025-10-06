<template>
  <!-- 移动端悬浮按钮 -->
  <button
    v-show="visible"
    class="toc-drawer-btn d-lg-none"
    @click="openDrawer"
    aria-label="打开本页目录"
    @touchstart="onTouchStart"
    @touchmove="onTouchMove"
    @touchend="onTouchEnd"
    :style="{ top: buttonTop + 'px' }"
  >
    <i class="fas fa-bookmark"></i>
  </button>
</template>

<script>
export default {
  name: 'TocDrawerButton',
  data() {
    return {
      visible: false,
      isDragging: false,
      startY: 0,
      initialTop: 0,
      buttonTop: window.innerHeight - 160,
      touchMoved: false
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
    window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: this.buttonTop + GAP } }));
  },
  beforeUnmount() {
    window.removeEventListener('scroll', this.onScroll);
    window.removeEventListener('resize', this.onResize);
    window.removeEventListener('floating-buttons-base-top', this.syncBaseTop);
  },
  methods: {
    // 接收“下方按钮”的基准位置，使 TOC 保持在其上方固定间距，且不越界
    syncBaseTop(e) {
      const base = e?.detail?.baseTop;
      const GAP = 48;
      const buttonHeight = 40;
      const minTop = 20;
      const maxTopUpper = window.innerHeight - buttonHeight - 20 - GAP; // 上方按钮的最大 top
      if (typeof base === 'number') {
        const desiredTop = base - GAP;
        this.buttonTop = Math.max(minTop, Math.min(maxTopUpper, desiredTop));
      }
    },
    onResize() {
      // 目录按钮在移动端常驻显示
      const isMobile = window.innerWidth < 992;
      this.visible = isMobile;
      // 初始位置靠近右下角，避开 BackToTop（默认更靠上）
      this.buttonTop = Math.min(window.innerHeight - 160, Math.max(20, this.buttonTop));
    },
    onScroll() {
      if (!this.isDragging) {
        const isMobile = window.innerWidth < 992;
        this.visible = isMobile;
      }
    },
    openDrawer() {
      if (!this.isDragging) {
        window.dispatchEvent(new Event('open-mobile-toc-drawer'));
      }
    },
    onTouchStart(e) {
      this.isDragging = true;
      this.touchMoved = false;
      this.startY = e.touches[0].clientY;
      this.initialTop = this.buttonTop;
    },
    onTouchMove(e) {
      this.touchMoved = true;
      const currentY = e.touches[0].clientY;
      const diffY = currentY - this.startY;
      const buttonHeight = 40;
      const GAP = 48;
      const minTop = 20;
      const maxTopUpper = window.innerHeight - buttonHeight - 20 - GAP; // 不把下方按钮推到底边
      let newTop = Math.max(minTop, Math.min(maxTopUpper, this.initialTop + diffY));
      this.buttonTop = newTop;

      // TOC 为上方按钮：派发“下方按钮”的基准位置，让 BackToTop 跟随在 GAP 下方
      window.dispatchEvent(new CustomEvent('floating-buttons-base-top', { detail: { baseTop: newTop + GAP } }));

      e.preventDefault();
    },
    onTouchEnd() {
      this.isDragging = false;
      if (!this.touchMoved) this.openDrawer();
    }
  }
};
</script>

<style scoped>
.toc-drawer-btn {
  position: fixed;
  right: 30px;
  width: 40px;
  height: 40px;
  background-color: #fff; /* 纯白 */
  color: var(--bs-body-color);
  border: 1px solid var(--bs-border-color);
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 18px;
  font-weight: bold;
  cursor: pointer;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12);
  z-index: 1000;
  outline: none;
  -webkit-tap-highlight-color: transparent;
  touch-action: none;
}
.toc-drawer-btn:hover {
  background-color: #fff; /* hover 仍保持纯白 */
  box-shadow: 0 4px 12px rgba(0, 0, 0, 0.18);
}
.toc-drawer-btn:active {
  transform: scale(0.95);
}
</style>