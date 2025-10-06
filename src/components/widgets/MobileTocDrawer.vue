<template>
  <!-- 仅移动端显示 -->
  <div v-if="show" class="mobile-offcanvas d-lg-none" @click.self="close">
    <div class="offcanvas-panel offcanvas-right border-start rounded-0 shadow-sm">
      <div class="offcanvas-header d-flex align-items-center justify-content-end">
        <button class="btn btn-sm p-2 btn-icon" @click="close" @focus="$event.target.blur()" aria-label="关闭">
          <i class="bi bi-x-lg"></i>
        </button>
      </div>

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
  name: 'MobileTocDrawer',
  components: { OnThisPage },
  data() {
    return {
      show: false,
      _openHandler: null
    };
  },
  mounted() {
    // 监听按钮触发事件
    this._openHandler = () => { this.show = true; };
    window.addEventListener('open-mobile-toc-drawer', this._openHandler);
  },
  beforeUnmount() {
    if (this._openHandler) {
      window.removeEventListener('open-mobile-toc-drawer', this._openHandler);
      this._openHandler = null;
    }
  },
  methods: {
    close() {
      this.show = false;
    }
  }
};
</script>

<style scoped>
/* 统一移动端抽屉样式（与 AppHeader 一致），右侧版本 */
.mobile-offcanvas {
  position: fixed;
  inset: 0;
  z-index: 1050;
}
.offcanvas-backdrop {
  position: absolute;
  inset: 0;
  background: rgba(0,0,0,0.25);
  z-index: 1;
}
.offcanvas-panel {
  position: absolute;
  top: 0; bottom: 0;
  width: min(85vw, 320px);
  background: #e9ecef;
  box-shadow: 0 0.25rem 0.75rem rgba(0,0,0,0.08);
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
  color: var(--bs-gray-700);
  margin-bottom: 0.25rem;
  font-weight: 600;
}
.offcanvas-card {
  background-color: #fff;
  border: none;
  border-radius: 0.5rem;
  padding: 0.5rem;
}
</style>